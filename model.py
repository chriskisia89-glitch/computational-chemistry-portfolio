"""Linear-in-parameters virial regression and thermodynamic calculations.

Canonical numerical references: N0=1 mmol/g, P0=1 pressure unit, T in K.
z=ln((P/P0)/(N/N0)) = sum(a_i*u**i)/T + sum(b_j*u**j).
a_i have kelvin units in this dimensionless-loading convention; b_j are
dimensionless. Physical predictions inherit the input pressure reference.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import brentq

R_KJ_MOL_K = 0.00831446261815324
PRESSURE_TO_TORR = {'Torr': 1.0, 'Pa': 760.0 / 101325.0,
                    'kPa': 760000.0 / 101325.0,
                    'bar': 76000000.0 / 101325.0,
                    'mbar': 76000.0 / 101325.0}
LOADING_TO_MMOL_G = {'mmol/g': 1.0, 'mol/kg': 1.0, 'mol/g': 1000.0}


def _vector(value, name, *, positive=False, nonnegative=False):
    out = np.asarray(value, dtype=float)
    if out.ndim != 1 or out.size == 0:
        raise ValueError(f'{name} must be a nonempty one-dimensional array')
    if not np.all(np.isfinite(out)):
        raise ValueError(f'{name} must contain only finite values')
    if positive and np.any(out <= 0):
        raise ValueError(f'{name} must be strictly positive')
    if nonnegative and np.any(out < 0):
        raise ValueError(f'{name} must be nonnegative')
    return out


def _order(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or not 0 <= value <= 8:
        raise ValueError(f'{name} must be an integer between 0 and 8')
    return int(value)


def design_matrix(loading, temperature, m=3, n=1):
    m, n = _order(m, 'm'), _order(n, 'n')
    u = _vector(loading, 'loading', positive=True)
    t = _vector(temperature, 'temperature', positive=True)
    if len(u) != len(t):
        raise ValueError('loading and temperature lengths differ')
    with np.errstate(over='raise', invalid='raise'):
        try:
            x = np.column_stack([u**i/t for i in range(m+1)] + [u**j for j in range(n+1)])
        except FloatingPointError as exc:
            raise ValueError('model matrix overflow') from exc
    if not np.all(np.isfinite(x)):
        raise ValueError('model matrix is not finite')
    return x


@dataclass(frozen=True)
class VirialFit:
    a: np.ndarray
    b: np.ndarray
    rank: int
    rmse_log: float
    r2_log: float | None
    condition_raw: float
    condition_scaled: float

    def __post_init__(self):
        for name in ('a', 'b'):
            v = _vector(getattr(self, name), name).copy()
            v.setflags(write=False)
            object.__setattr__(self, name, v)

    def log_p_over_n(self, loading, temperature):
        u, t = np.broadcast_arrays(np.asarray(loading, dtype=float), np.asarray(temperature, dtype=float))
        if not np.all(np.isfinite(u)) or np.any(u < 0):
            raise ValueError('prediction loading must be finite and nonnegative')
        if not np.all(np.isfinite(t)) or np.any(t <= 0):
            raise ValueError('prediction temperature must be finite and positive')
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            try:
                z = np.polynomial.polynomial.polyval(u, self.a)/t + np.polynomial.polynomial.polyval(u, self.b)
            except FloatingPointError as exc:
                raise ValueError('prediction overflow') from exc
        if not np.all(np.isfinite(z)):
            raise ValueError('nonfinite prediction')
        return z

    def pressure(self, loading, temperature):
        u, t = np.broadcast_arrays(np.asarray(loading, dtype=float), np.asarray(temperature, dtype=float))
        z = self.log_p_over_n(u, t)
        out = np.zeros(u.shape, dtype=float)
        positive = u > 0
        try:
            with np.errstate(over='raise', invalid='raise', under='ignore'):
                out[positive] = np.exp(np.log(u[positive]) + z[positive])
        except FloatingPointError as exc:
            raise ValueError('pressure overflow') from exc
        if np.any(positive & (out == 0)):
            raise ValueError('pressure underflow')
        return out

    def heat(self, loading):
        u = np.asarray(loading, dtype=float)
        if not np.all(np.isfinite(u)) or np.any(u < 0):
            raise ValueError('heat loading must be finite and nonnegative')
        with np.errstate(over='raise', invalid='raise'):
            try:
                q = -R_KJ_MOL_K * np.polynomial.polynomial.polyval(u, self.a)
            except FloatingPointError as exc:
                raise ValueError('heat calculation overflow') from exc
        if not np.all(np.isfinite(q)):
            raise ValueError('nonfinite heat')
        return q

    @property
    def q0(self):
        return float(-R_KJ_MOL_K*self.a[0])

    def inverse_loading(self, pressure, temperature, bracket):
        p, t = float(pressure), float(temperature)
        lo, hi = map(float, bracket)
        if not np.isfinite([p, t, lo, hi]).all() or p <= 0 or t <= 0 or lo < 0 or hi <= lo:
            raise ValueError('positive finite pressure/temperature and 0 <= lo < hi are required')
        # sign(dP/du) = sign(1 + u*(A'(u)/T + B'(u))) for u>0.
        # Check the polynomial at its endpoints AND every internal stationary point.
        g = Polynomial([1.0]) + Polynomial([0.0, 1.0])*(Polynomial(self.a).deriv()/t + Polynomial(self.b).deriv())
        locations = [lo, hi]
        for root in g.deriv().roots():
            if abs(root.imag) < 1e-10 and lo < root.real < hi:
                locations.append(float(root.real))
        if np.min(g(locations)) <= 0:
            raise ValueError('pressure is not strictly increasing throughout the requested bracket')
        p_lo, p_hi = float(self.pressure(lo, t)), float(self.pressure(hi, t))
        if not p_lo <= p <= p_hi:
            raise ValueError('pressure lies outside the requested loading bracket')
        return float(brentq(lambda u: float(self.pressure(u, t))-p, lo, hi, xtol=1e-12, rtol=1e-12))


def fit_virial(loading, temperature, log_p_over_n, m=3, n=1):
    x = design_matrix(loading, temperature, m, n)
    z = _vector(log_p_over_n, 'log_p_over_n')
    if len(z) != len(x):
        raise ValueError('response and model-matrix lengths differ')
    if x.shape[0] < x.shape[1]:
        raise ValueError('insufficient observations for the requested model')
    scales = np.linalg.norm(x, axis=0)
    if np.any(scales == 0) or not np.all(np.isfinite(scales)):
        raise ValueError('model has invalid column scales')
    scaled = x/scales
    beta_scaled, _, rank, _ = np.linalg.lstsq(scaled, z, rcond=None)
    if rank < x.shape[1]:
        raise ValueError('rank-deficient fit: use distinct temperatures and adequate loading variation')
    beta = beta_scaled/scales
    residual = z-x@beta
    sse = float(residual@residual)
    sst = float(np.sum((z-z.mean())**2))
    return VirialFit(beta[:m+1], beta[m+1:], int(rank), float(np.sqrt(sse/len(z))),
                     None if sst == 0 else float(1-sse/sst),
                     float(np.linalg.cond(x)), float(np.linalg.cond(scaled)))


def load_pressure_csv(path: str | Path, *, pressure_unit='Torr', loading_unit='mmol/g'):
    if pressure_unit not in PRESSURE_TO_TORR or loading_unit not in LOADING_TO_MMOL_G:
        raise ValueError('unsupported pressure or loading unit')
    with Path(path).open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    try:
        temperature = _vector([float(r['temperature_K']) for r in rows], 'temperature', positive=True)
        loading = _vector([float(r['uptake']) for r in rows], 'uptake', positive=True)*LOADING_TO_MMOL_G[loading_unit]
        pressure = _vector([float(r['pressure']) for r in rows], 'pressure', positive=True)*PRESSURE_TO_TORR[pressure_unit]
    except (KeyError, TypeError) as exc:
        raise ValueError('CSV requires temperature_K, uptake and pressure columns') from exc
    return loading, temperature, np.log(pressure)-np.log(loading)


def benchmark_answers(path):
    u, t, z = load_pressure_csv(path)
    fit = fit_virial(u, t, z, 3, 1)
    loadings = np.array([0.1, 0.3, 0.5])
    return {'coefficients': {'a': fit.a.tolist(), 'b': fit.b.tolist()},
            'q0_kj_mol': fit.q0,
            'qst_kj_mol_at_loadings': fit.heat(loadings).tolist(),
            'pressure_torr': {str(temp): fit.pressure(loadings, temp).tolist() for temp in (273, 298)},
            'uptake_mmol_g_at_273K_for_pressure_torr': [fit.inverse_loading(p, 273, (float(u[t==273].min()), float(u[t==273].max()))) for p in (10, 50, 100)],
            'rmse_log': fit.rmse_log, 'rank': fit.rank}

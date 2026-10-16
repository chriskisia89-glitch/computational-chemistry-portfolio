# Two-temperature virial adsorption calculation

Implement `solve(csv_path)` in `submission_template.py`. Return a Python
dictionary with the numerical outputs below. The supplied CSV is an explicitly
synthetic dataset for a scientific coding benchmark, not experimental data.

## Input and model

The CSV columns are `temperature_K`, `uptake` in mmol/g and `pressure` in Torr.
Rows occur at 273 and 298 K and may appear in any order. All input numbers must
be finite and strictly positive.

Let `u = N / (1 mmol/g)` and `p = P / (1 Torr)`. Use natural logarithms and
fit all rows simultaneously to

```text
z = ln(p/u)
z = (a0 + a1*u + a2*u^2 + a3*u^3)/T + b0 + b1*u.
```

Estimate the six coefficients by **unweighted least squares in z**. The a
coefficients have kelvin units in this normalized-loading convention; the b
coefficients are dimensionless. Use the column order
`[1/T, u/T, u^2/T, u^3/T, 1, u]`. Do not fit residuals in pressure instead.
Use an identifiable full-rank system. SVD or QR with suitable column scaling
is recommended; an underdetermined or rank-deficient fit must raise ValueError.

Take `R = 0.00831446261815324 kJ/(mol K)`. Compute

```text
Qst(u) = -R*(a0 + a1*u + a2*u^2 + a3*u^3)
Q0 = -R*a0
P(u,T) = u*exp((a0+a1*u+a2*u^2+a3*u^3)/T+b0+b1*u) Torr.
```

Q0 is a zero-coverage mathematical limit. It is not an observed data point.

## Required return values

| Key | Value |
| --- | --- |
| `coefficients` | `{"a": [a0,a1,a2,a3], "b": [b0,b1]}` |
| `q0_kj_mol` | Q0 in kJ/mol |
| `qst_kj_mol_at_loadings` | Qst at N = [0.1, 0.3, 0.5] mmol/g, in that order |
| `pressure_torr` | `{"273": [P(0.1,273), P(0.3,273), P(0.5,273)], "298": [P(0.1,298), P(0.3,298), P(0.5,298)]}` |
| `uptake_mmol_g_at_273K_for_pressure_torr` | Uptakes at P = [10,50,100] Torr, T = 273 K, in that order |
| `rmse_log` | sqrt(mean((z_observed-z_fitted)^2)) over all rows |
| `rank` | Integer rank of the fitted six-column system |

Find each inverse uptake within the minimum-to-maximum **273 K loading range**
in the input. Check that the model is strictly increasing throughout that
bracket and that the target pressure is bracketed; otherwise raise ValueError.
Do not extrapolate an inverse outside the bracket. The two evaluation datasets
have valid monotonic inverses for the requested targets.

## Input validation

Raise ValueError for missing required columns, nonfinite values, nonpositive
pressure, uptake or temperature, fewer than six usable rows, and a
rank-deficient design. A single-temperature dataset cannot uniquely identify
the shared temperature-dependent coefficients in this specification.

## Numerical acceptance tolerances

These apply to every scalar in the corresponding output, with the usual
`abs(actual-expected) <= max(atol, rtol*abs(expected))` interpretation.

| Output | Relative tolerance | Absolute tolerance |
| --- | --- | --- |
| Coefficients | 1e-7 | 1e-5 |
| Q0 and Qst | 1e-7 | 1e-7 kJ/mol |
| Pressures | 1e-7 | 1e-5 Torr |
| Inverse uptakes | 1e-7 | 1e-6 mmol/g |
| Log RMSE | 1e-6 | 1e-9 |
| Rank | exact integer | exact integer |

All returned numbers must be finite. Preserve the exact keys, list lengths,
and orders above. Do not include additional top-level keys.

## Scoring

Two numerical datasets are evaluated, and each contributes half of the
following component weights: coefficients 25, heats 20, pressures 20,
inverse uptakes 15, and RMSE/rank 10. Eight invalid-input cases contribute
1.25 points each. Total: 100. The second numerical dataset and frozen answers
are retained in the instructor package, outside the candidate-only archive.

The full instructor package contains a reference implementation, evaluator,
frozen QR-based numerical answers and implementation tests. The candidate
archive contains this specification, the supplied CSV, template and dependency
list; it omits solutions, generation truth and evaluator fixtures.

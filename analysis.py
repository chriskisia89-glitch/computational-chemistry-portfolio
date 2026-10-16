"""Run the published-plot audit and controlled synthetic analysis."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.linalg import lstsq as qr_lstsq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .model import R_KJ_MOL_K, design_matrix, fit_virial, load_pressure_csv, benchmark_answers


def _save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def audit(root):
    root=Path(root);out=root/'results';out.mkdir(exist_ok=True)
    with (root/'data/published/c_pbdaz_300_digitized.csv').open() as f:
        rows=list(csv.DictReader(f))
    u=np.array([float(r['uptake_mmol_g']) for r in rows])
    t=np.array([float(r['temperature_K']) for r in rows])
    z=np.array([float(r['log_p_over_n_source']) for r in rows])
    fit=fit_virial(u,t,z)
    published=json.loads((root/'data/published/published_coefficients.json').read_text())
    pub_a=np.array(published['a']);pub_b=np.array(published['b'])
    x=design_matrix(u,t)
    pub_z=x@np.r_[pub_a,pub_b]
    fitted_z=fit.log_p_over_n(u,t)
    offset=float(np.mean(z-pub_z))
    pub_rmse=float(np.sqrt(np.mean((z-pub_z)**2)))
    shifted_rmse=float(np.sqrt(np.mean((z-pub_z-offset)**2)))
    # Independent QR algorithm check, compared with the scaled NumPy SVD fit.
    scale=np.linalg.norm(x,axis=0)
    qr_beta=qr_lstsq(x/scale,z,lapack_driver='gelsy')[0]/scale
    qr_delta=float(np.max(np.abs(qr_beta-np.r_[fit.a,fit.b])))
    common_lo=max(u[t==273].min(),u[t==298].min())
    common_hi=min(u[t==273].max(),u[t==298].max())
    report={
      'data_type':'approximate digitized coordinates from a published plot; not raw measurements',
      'rows':len(rows),'counts_by_temperature':{str(int(temp)):int(np.sum(t==temp)) for temp in np.unique(t)},
      'reference_equation':'ln((P/P0)/(N/N0)) = sum(a_i*(N/N0)^i)/T + sum(b_j*(N/N0)^j); m=3,n=1',
      'source_pressure_reference':'unspecified in the selected panel; no Torr/bar label assigned to inferred source pressures',
      'fit_a_K':fit.a.tolist(),'fit_b':fit.b.tolist(),'fit_rmse_log':fit.rmse_log,
      'fit_r2_log':fit.r2_log,'rank':fit.rank,
      'condition_number_raw':fit.condition_raw,'condition_number_column_scaled':fit.condition_scaled,
      'independent_qr_max_coefficient_delta':qr_delta,
      'published_a_K':pub_a.tolist(),'published_b':pub_b.tolist(),
      'published_coefficients_rmse_against_digitized_logordinate':pub_rmse,
      'best_global_log_offset':offset,
      'published_coefficients_rmse_after_best_global_log_offset':shifted_rmse,
      'published_a0_implied_q0_kj_mol':float(-R_KJ_MOL_K*pub_a[0]),
      'published_rounded_q0_kj_mol':published['reported_q0_kj_mol_rounded'],
      'digitized_refit_extrapolated_q0_kj_mol':fit.q0,
      'common_observed_loading_interval_mmol_g':[float(common_lo),float(common_hi)],
      'interpretation':'The printed coefficients and digitized panel disagree under the printed equation, even after one global log offset. This project does not resolve the cause, reproduce the authors\' original fitting procedure, or treat the graphical refit heat as an independently validated physical measurement.',
      'zero_coverage_note':'Q0 is an extrapolation to N=0; it is not an observed data point.',
      'source_license':'CC BY-NC 4.0; see THIRD_PARTY_NOTICES.md'}
    _save_json(out/'published_audit.json',report)
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors={273:'#2455a4',298:'#bb3c36'}
    for temp in (273,298):
        chosen=t==temp;grid=np.linspace(u[chosen].min(),u[chosen].max(),250)
        axs[0,0].scatter(u[chosen],z[chosen],s=20,color=colors[temp],label=f'{temp} K digitized')
        axs[0,0].plot(grid,fit.log_p_over_n(grid,temp),color=colors[temp],label=f'{temp} K refit')
        pub_grid=np.polynomial.polynomial.polyval(grid,pub_a)/temp+np.polynomial.polynomial.polyval(grid,pub_b)
        axs[0,0].plot(grid,pub_grid,'--',color=colors[temp],alpha=.65,label=f'{temp} K printed coefficients')
        axs[0,1].scatter(u[chosen],z[chosen]-fitted_z[chosen],s=20,color=colors[temp],label=f'{temp} K')
        axs[1,1].scatter(u[chosen],z[chosen]-pub_z[chosen]-offset,s=20,color=colors[temp],label=f'{temp} K')
    axs[0,0].set(title='Published plot: digitization, refit and printed coefficients',xlabel='Uptake (mmol/g)',ylabel='Source ln(P/N) ordinate')
    axs[0,0].legend(fontsize=7,ncol=2)
    axs[0,1].axhline(0,color='black',lw=.6);axs[0,1].set(title='Residuals of graphical refit',xlabel='Uptake (mmol/g)',ylabel='Digitized minus refit (log ordinate)')
    heat_grid=np.linspace(common_lo,common_hi,150)
    axs[1,0].plot(heat_grid,fit.heat(heat_grid),label='Graphical refit: mathematical result',color='#1b7a57')
    axs[1,0].plot(heat_grid,-R_KJ_MOL_K*np.polynomial.polynomial.polyval(heat_grid,pub_a),label='Printed coefficient calculation',color='#7c5295')
    axs[1,0].set(title='Derived heat over shared observed loading range',xlabel='Uptake (mmol/g)',ylabel='Qst (kJ/mol)');axs[1,0].legend(fontsize=7)
    axs[1,1].axhline(0,color='black',lw=.6);axs[1,1].set(title='Printed-coefficient discrepancy after global offset',xlabel='Uptake (mmol/g)',ylabel='Remaining log-ordinate discrepancy')
    for ax in axs.flat:ax.grid(alpha=.2)
    fig.suptitle('C-PBDAZ-300 reproducibility audit | digitized source data',fontsize=13)
    fig.savefig(out/'published_audit.png',dpi=180);plt.close(fig)
    with (out/'published_fit_predictions.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['temperature_K','uptake_mmol_g','digitized_logordinate','refit_logordinate','printed_coefficient_logordinate'])
        w.writerows(zip(t,u,z,fitted_z,pub_z))
    return report


def controlled(root):
    root=Path(root);out=root/'results';out.mkdir(exist_ok=True)
    answers=benchmark_answers(root/'task/isotherms.csv')
    _save_json(out/'benchmark_answers.json',answers)
    u,t,z=load_pressure_csv(root/'task/isotherms.csv')
    fit=fit_virial(u,t,z)
    x=design_matrix(u,t);sc=np.linalg.norm(x,axis=0)
    qr=qr_lstsq(x/sc,z,lapack_driver='gelsy')[0]/sc
    _save_json(out/'numerical_crosscheck.json',{
        'numpy_method':'column-scaled SVD least squares',
        'scipy_method':'column-scaled QR least squares with gelsy',
        'max_coefficient_absolute_difference':float(np.max(np.abs(qr-np.r_[fit.a,fit.b]))),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__})
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for temp,color in [(273,'#2455a4'),(298,'#bb3c36')]:
        selected=t==temp;grid=np.linspace(u[selected].min(),u[selected].max(),200)
        axs[0].scatter(u[selected],u[selected]*np.exp(z[selected]),s=17,color=color,label=f'{temp} K synthetic data')
        axs[0].plot(grid,fit.pressure(grid,temp),color=color,label=f'{temp} K fit')
        axs[1].scatter(u[selected],z[selected]-fit.log_p_over_n(u[selected],temp),s=18,color=color,label=f'{temp} K')
    axs[0].set(xlabel='Uptake (mmol/g)',ylabel='Pressure (Torr)',title='Synthetic benchmark: fitted isotherms')
    axs[1].set(xlabel='Uptake (mmol/g)',ylabel='Residual in ln(P/N)',title='Synthetic benchmark: fitting residuals');axs[1].axhline(0,color='black',lw=.6)
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(out/'synthetic_benchmark.png',dpi=180);plt.close(fig)
    return answers


def write_report(root,published,benchmark):
    root=Path(root)
    md=f'''# Adsorption analysis: results and limits

Created on 16 October 2026 as an AI-assisted portfolio project.

## 1. Controlled coding benchmark

The input contains 70 **synthetic** observations at 273 and 298 K. It is
generated by `scripts/generate_benchmark.py` with a fixed random seed and is
not a reproduction of the original experiments. The fitting objective is
unweighted least squares in the transformed log ordinate, not in pressure.

- Matrix rank: {benchmark['rank']} of 6.
- Log-ordinate RMSE: {benchmark['rmse_log']:.10g}.
- Fitted zero-coverage heat: {benchmark['q0_kj_mol']:.10g} kJ/mol.
- Coefficients and all requested predictions: `benchmark_answers.json`.
- Independent SVD/QR comparison: `numerical_crosscheck.json`.

The heat is calculated as Qst(u) = -R sum(a_i u^i), with
R = {R_KJ_MOL_K:.17g} kJ/(mol K). The source equation uses numerical
loading normalized by 1 mmol/g and temperature in kelvin. Q0 = -R a0 is
an extrapolated limiting value, not a measured zero-loading observation.

## 2. Published-plot audit

Source: Edlabadkar et al., Chemistry of Materials (2024),
DOI https://doi.org/10.1021/acs.chemmater.3c01717;
Supporting Information DOI https://doi.org/10.1021/acs.chemmater.3c01717.s001,
S-5, Table S.2, C-PBDAZ-300.

The {published['rows']} extracted coordinates comprise
{published['counts_by_temperature']['273']} at 273 K and
{published['counts_by_temperature']['298']} at 298 K. They are approximate
coordinates of plotted markers/curves, not the original measurement file.
Dense red markers are thinned by an explicit image-processing rule.
The source image is audited by checksum and axis-tick calibration in
`data/published/provenance.json`; `digitization_audit.png` shows the selection.

- Graphical-refit log-ordinate RMSE: {published['fit_rmse_log']:.10g}.
- RMSE using the printed coefficients: {published['published_coefficients_rmse_against_digitized_logordinate']:.10g}.
- RMSE after the best common log offset: {published['published_coefficients_rmse_after_best_global_log_offset']:.10g}.
- Q0 calculated from the printed a0: {published['published_a0_implied_q0_kj_mol']:.10g} kJ/mol;
  the source reports the rounded value 21 kJ/mol.
- Q0 from the graphical refit: {published['digitized_refit_extrapolated_q0_kj_mol']:.10g} kJ/mol.

The printed coefficients do not reproduce the selected raster plot under its
printed equation. One constant log offset, such as a change of pressure
reference, does not eliminate the discrepancy. Its cause cannot be established
from this panel. The graphical refit does not establish corrected physical heat
values and is not an exact reproduction of the original fitting procedure,
published R2, reduced chi-square or experimental uncertainties.

The source panel does not state the pressure reference unit. Its displayed
log ordinate is therefore fitted directly. No absolute source pressures are
labelled as Torr or bar. Only the independent synthetic task explicitly uses Torr.
Heat curves are displayed over the shared observed loading interval
{published['common_observed_loading_interval_mmol_g'][0]:.6g} to
{published['common_observed_loading_interval_mmol_g'][1]:.6g} mmol/g.

## 3. What has and has not been validated

The automated tests assess known coefficient recovery, independent analytic
cases, thermodynamic derivatives, unit conversion, rank deficiency, invalid
inputs, inversion and grading behavior. The plot refit is descriptive; pixel
resolution, merged markers and fit lines prevent experimental uncertainty
estimation. Numerical agreement on synthetic data validates implementation,
not the physical adequacy of a model for every adsorbent.

## 4. Attribution and reuse

Original experiments belong to the cited researchers. This new code and
analysis were prepared with OpenAI assistance at the user's request; the
package does not establish the user's historical contribution to those studies.
The source-derived data and audit images retain the source CC BY-NC 4.0
designation. The new code and independent synthetic task are under the MIT
license. See `THIRD_PARTY_NOTICES.md` for the division of materials.
'''
    (root/'results/REPORT.md').write_text(md,encoding='utf-8')


def run(root):
    root=Path(root);p=audit(root);b=controlled(root);write_report(root,p,b)
    print(json.dumps({'source_rows':p['rows'],'synthetic_rows':70,
      'synthetic_rmse_log':b['rmse_log'],'synthetic_q0_kj_mol':b['q0_kj_mol'],
      'source_printed_q0_kj_mol':p['published_a0_implied_q0_kj_mol']},indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    run(parser.parse_args().root)

# Virial adsorption analysis and coding-task portfolio

A runnable Python project created on 16 October 2026 with OpenAI assistance.
It combines an audit of published CO2 adsorption plots with a separate,
fully controlled scientific coding task. Original experiments remain credited
to the published researchers. This package is new portfolio work and does not
verify anyone's past research contribution, degree or employment history.

## Start here

Use Python 3.11 or newer; this release was tested with Python 3.12.

### Windows

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python scripts\run_all.py
```

### macOS or Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/run_all.py
```

The normal run is offline after dependency installation. It uses the included
CSV data, recreates plots and results, evaluates the reference implementation,
and runs pytest. It does not rewrite the frozen reference answers or regenerate
input data. A failing calculation, evaluation or test produces a nonzero exit.

## What is included

| Location | Contents |
| --- | --- |
| `src/virial_portfolio/model.py` | Column-scaled regression, heat calculation, pressure prediction, inversion and unit handling |
| `data/published/` | Digitized source-plot coordinates, provenance and printed coefficients |
| `task/isotherms.csv` | 70 explicitly synthetic observations at 273 and 298 K |
| `task/TASK.md` | Standalone coding question, equations, outputs, validation rules and tolerances |
| `task/reference_solution.py` | Runnable reference implementation |
| `task/expected_answers.json` | Frozen answers from an independent QR calculation |
| `task/evaluate_solution.py` | Numerical and invalid-input evaluator with component scores |
| `task/grade_submission.py` | JSON-only numerical comparison utility |
| `tests/` | Analytic cases, unit checks, regression checks and grader tests |
| `results/` | Report, figures, numerical answers and actual test/evaluation logs |
| `candidate_bundle.zip` | Candidate-facing question, data and template, without solutions |
| `APPLICATION_NOTES.md` | Accurate description of this new AI-assisted portfolio work |

## Two distinct kinds of evidence

The published-data portion uses approximate coordinates from the raster plot
for C-PBDAZ-300 in Table S.2, Supporting Information S-5, Edlabadkar et al.
(2024). These are not raw instrument measurements. The selected plot and its
printed coefficients disagree under the printed model. The report preserves
that observation and does not present the graphical refit as corrected
physical measurements or an exact reproduction of the original authors' fit.
The selected panel does not state the pressure reference unit, so source log
ordinates are fitted directly without assigning Torr or bar to source pressures.

The synthetic task uses independently chosen parameters, specified Torr and
mmol/g units, and controlled noise. It has frozen numerical answers, a second
validation dataset and a reference evaluator. It tests computation without
depending on access to a paper or on the graphical discrepancy.

Read `results/REPORT.md` and `results/validation_summary.json` for results,
limits and the number of tests actually run. A numerical implementation can
pass tests even when a published plot cannot be reproduced; these are separate
scientific conclusions.

## Run a reference or candidate solution

```bash
python task/reference_solution.py task/isotherms.csv results/reference_answers.json
python task/grade_submission.py results/reference_answers.json
python task/evaluate_solution.py task/reference_solution.py --output results/reference_evaluation.json
```

For another implementation, replace `task/reference_solution.py` in the final
command with the path to its Python module containing `solve(csv_path)`.
The evaluator imports and executes that module. The JSON-only grading command
does not execute candidate code.

## Reproduce the source digitization

This optional step needs `requirements-source.txt`. Obtain the Supporting
Information via DOI https://doi.org/10.1021/acs.chemmater.3c01717.s001 or
the file link in `data/published/provenance.json`, then run:

```bash
python scripts/extract_published_data.py /path/to/cm3c01717_si_001.pdf
```

The extractor checks the PDF checksum, extracts its original raster, calibrates
axis ticks, detects colored marker/curve centers and writes the selected
coordinates and audit image. It fails rather than silently processing a
different file version. The PDF itself is not redistributed in this package.

`scripts/generate_benchmark.py` and `scripts/freeze_reference_answers.py` are
maintainer tools. Ordinary validation does not run them. Changes to inputs,
model choices or golden answers require review and a new results run.

## Attribution and distribution

The new code and standalone synthetic task use the MIT license in `LICENSE`.
Source-derived coordinates and audit images retain the supporting record's
CC BY-NC 4.0 designation; see `THIRD_PARTY_NOTICES.md`. The candidate-only
archive contains no source-derived data or images. Keep the attribution and
license notices with the corresponding materials.

The project has not been published to GitHub. To obtain a repository URL,
upload the project to your own account after reviewing its contents and the
component licenses. Describe your personal contribution accurately.

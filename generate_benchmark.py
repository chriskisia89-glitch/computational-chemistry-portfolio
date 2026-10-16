"""Generate an explicitly synthetic task dataset, independent of paper data."""
from pathlib import Path
import csv
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def generate():
    rng=np.random.default_rng(271828)
    a=np.array([-2500.,500.,-30.,1.5]);b=np.array([12.,-1.])
    rows=[]
    for temp,lo,hi,count in [(273,.05,2.5,40),(298,.02,.7,30)]:
        for u in np.linspace(lo,hi,count):
            # Independent explicit formula, not a call to the model under test.
            z=(-2500+500*u-30*u*u+1.5*u*u*u)/temp+12-u
            noise=float(rng.normal(0,.004))
            p=u*np.exp(z+noise)
            rows.append({'temperature_K':temp,'uptake':f'{u:.16g}','pressure':f'{p:.16g}'})
    dest=ROOT/'task';dest.mkdir(exist_ok=True)
    with (dest/'isotherms.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['temperature_K','uptake','pressure']);w.writeheader();w.writerows(rows)
    (dest/'generation_truth.json').write_text(json.dumps({'data_type':'synthetic, not experimental measurements',
      'seed':271828,'noise':'Gaussian additive noise in ln(P/N), standard deviation 0.004',
      'pressure_unit':'Torr','loading_unit':'mmol/g','temperatures_K':[273,298],
      'true_a':a.tolist(),'true_b':b.tolist(),'rows':len(rows),
      'purpose':'Controlled coding benchmark with known generating coefficients; no data borrowed from the paper.'},indent=2)+'\n')
    print(f'Generated {len(rows)} synthetic rows')
    # A second dataset is reserved for evaluator validation and omitted from
    # the candidate-only archive. It prevents scoring a one-file answer lookup
    # as a generally working calculation.
    variant=[]
    for temp,lo,hi,count in [(273,.05,2.5,31),(298,.02,.7,23)]:
        for u in np.linspace(lo,hi,count):
            z=(-2400+550*u-35*u*u+u**3)/temp+13-1.1*u
            variant.append([temp,f'{u:.16g}',f'{u*np.exp(z):.16g}'])
    fixtures=ROOT/'tests/fixtures';fixtures.mkdir(parents=True,exist_ok=True)
    with (fixtures/'validation_isotherms.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['temperature_K','uptake','pressure']);w.writerows(variant)

if __name__=='__main__':generate()

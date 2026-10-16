"""Maintainer action: compute frozen answers with independent QR/scalar formulas.

Not part of the ordinary run/test pipeline. Updating golden answers requires
deliberate scientific review; tests must not regenerate their own expectations.
"""
from pathlib import Path
import csv,json,math
import numpy as np
from scipy.linalg import lstsq
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]


def independent_answers(path):
    with Path(path).open() as f:r=list(csv.DictReader(f))
    u=np.array([float(x['uptake']) for x in r]);t=np.array([float(x['temperature_K']) for x in r])
    z=np.array([math.log(float(x['pressure']))-math.log(float(x['uptake'])) for x in r])
    x=np.array([[v**0/temp,v/temp,v*v/temp,v**3/temp,1,v] for v,temp in zip(u,t)])
    scales=np.sqrt(np.sum(x*x,axis=0));c=lstsq(x/scales,z,lapack_driver='gelsy')[0]/scales
    a0,a1,a2,a3,b0,b1=map(float,c)
    def p(v,temp):return v*math.exp((a0+a1*v+a2*v*v+a3*v*v*v)/temp+b0+b1*v)
    def q(v):return -.00831446261815324*(a0+a1*v+a2*v*v+a3*v*v*v)
    return {'coefficients':{'a':[a0,a1,a2,a3],'b':[b0,b1]},
      'q0_kj_mol':q(0),'qst_kj_mol_at_loadings':[q(v) for v in [.1,.3,.5]],
      'pressure_torr':{str(temp):[p(v,temp) for v in [.1,.3,.5]] for temp in [273,298]},
      'uptake_mmol_g_at_273K_for_pressure_torr':[float(brentq(lambda v:p(v,273)-target,float(u[t==273].min()),float(u[t==273].max()),xtol=1e-13)) for target in [10,50,100]],
      'rmse_log':float(np.sqrt(np.mean((z-x@c)**2))),'rank':int(np.linalg.matrix_rank(x/scales))}


if __name__=='__main__':
    for source,dest in [(ROOT/'task/isotherms.csv',ROOT/'task/expected_answers.json'),(ROOT/'tests/fixtures/validation_isotherms.csv',ROOT/'tests/fixtures/validation_expected.json')]:
        dest.write_text(json.dumps(independent_answers(source),indent=2,allow_nan=False)+'\n')
        print(dest)

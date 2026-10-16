"""Compare a JSON answer with a frozen numerical reference; no candidate code execution."""
import argparse
import json
import math
from pathlib import Path

TOLERANCES={
 'coefficients':(1e-7,1e-5), 'q0_kj_mol':(1e-7,1e-7),
 'qst_kj_mol_at_loadings':(1e-7,1e-7), 'pressure_torr':(1e-7,1e-5),
 'uptake_mmol_g_at_273K_for_pressure_torr':(1e-7,1e-6),
 'rmse_log':(1e-6,1e-9),'rank':(0,0)}


def compare(actual,expected,rtol,atol,path='answer'):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or set(actual)!=set(expected):return [f'{path}: incorrect keys']
        return sum((compare(actual[k],v,rtol,atol,f'{path}.{k}') for k,v in expected.items()),[])
    if isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):return [f'{path}: incorrect list length']
        return sum((compare(a,e,rtol,atol,f'{path}[{i}]') for i,(a,e) in enumerate(zip(actual,expected))),[])
    if isinstance(actual,bool) or not isinstance(actual,(int,float)) or not math.isfinite(actual):
        return [f'{path}: expected a finite number']
    if abs(actual-expected)>max(atol,rtol*abs(expected)):return [f'{path}: outside numerical tolerance']
    return []


def validate_answer(actual,expected):
    if not isinstance(actual,dict) or set(actual)!=set(expected):return ['answer: incorrect top-level keys']
    errors=[]
    for key,value in expected.items():
        errors.extend(compare(actual[key],value,*TOLERANCES[key],key))
    if not isinstance(actual.get('rank'),int) or isinstance(actual.get('rank'),bool):errors.append('rank: expected an integer')
    return errors


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('answer_json');p.add_argument('--expected',default=str(Path(__file__).with_name('expected_answers.json')))
    a=p.parse_args();errors=validate_answer(json.loads(Path(a.answer_json).read_text()),json.loads(Path(a.expected).read_text()))
    print(json.dumps({'passed':not errors,'errors':errors},indent=2));raise SystemExit(bool(errors))

"""Evaluate solve(csv_path) on two numerical datasets and invalid-input cases."""
from pathlib import Path
import argparse,csv,importlib.util,json,tempfile,sys
from grade_submission import validate_answer, compare, TOLERANCES

ROOT=Path(__file__).resolve().parents[1]


def evaluate(solution_path):
    spec=importlib.util.spec_from_file_location('candidate_solution',solution_path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    results=[]
    with tempfile.TemporaryDirectory() as folder:
        for label,source,gold in [('provided',ROOT/'task/isotherms.csv',ROOT/'task/expected_answers.json'),('validation',ROOT/'tests/fixtures/validation_isotherms.csv',ROOT/'tests/fixtures/validation_expected.json')]:
            copy=Path(folder)/f'{label}.csv';copy.write_bytes(source.read_bytes())
            groups={'coefficients':(['coefficients'],25),'heat':(['q0_kj_mol','qst_kj_mol_at_loadings'],20),'pressure':(['pressure_torr'],20),'inverse':(['uptake_mmol_g_at_273K_for_pressure_torr'],15),'fit_statistics':(['rmse_log','rank'],10)}
            group_points={name:0 for name in groups}
            try:
                expected=json.loads(gold.read_text());actual=module.solve(copy)
                errors=validate_answer(actual,expected)
                if isinstance(actual,dict) and set(actual)==set(expected):
                    for group,(keys,weight) in groups.items():
                        issues=sum((compare(actual[k],expected[k],*TOLERANCES[k],k) for k in keys),[])
                        if 'rank' in keys and (not isinstance(actual['rank'],int) or isinstance(actual['rank'],bool)):issues.append('rank must be an integer')
                        if not issues:group_points[group]=weight/2
            except Exception as exc:errors=[f'{type(exc).__name__}: {exc}']
            results.append({'case':label,'passed':not errors,'errors':errors,'points':sum(group_points.values()),'component_points':group_points})
        with (ROOT/'task/isotherms.csv').open() as f:base=list(csv.DictReader(f))
        cases={}
        for name,key,value in [('zero_pressure','pressure','0'),('negative_loading','uptake','-1'),('nonfinite_pressure','pressure','nan'),('zero_temperature','temperature_K','0')]:
            rows=[r.copy() for r in base];rows[0][key]=value;cases[name]=rows
        cases['single_temperature']=[dict(r,temperature_K='273') for r in base]
        cases['insufficient_rows']=base[:5]
        cases['constant_loading']=[dict(r,uptake='1') for r in base]
        cases['missing_column']=[{k:v for k,v in r.items() if k!='pressure'} for r in base]
        for name,rows in cases.items():
            p=Path(folder)/f'{name}.csv'
            with p.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            passed=False;error='expected ValueError'
            try:module.solve(p)
            except ValueError:passed=True;error=''
            except Exception as exc:error=f'expected ValueError, got {type(exc).__name__}'
            results.append({'case':name,'passed':passed,'errors':[] if passed else [error],'points':10/len(cases) if passed else 0})
    return {'score_out_of_100':sum(r['points'] for r in results),'passed_all_cases':all(r['passed'] for r in results),'cases':results}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('solution_py');p.add_argument('--output')
    a=p.parse_args();result=evaluate(a.solution_py);txt=json.dumps(result,indent=2,allow_nan=False)
    if a.output:Path(a.output).write_text(txt+'\n')
    print(txt);raise SystemExit(not result['passed_all_cases'])

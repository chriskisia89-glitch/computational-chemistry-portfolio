"""Recreate analysis outputs and run tests without changing frozen references."""
from pathlib import Path
import json,os,subprocess,sys,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from virial_portfolio.analysis import run


def main():
    run(ROOT)
    evaluation=subprocess.run([sys.executable,str(ROOT/'task/evaluate_solution.py'),str(ROOT/'task/reference_solution.py'),'--output',str(ROOT/'results/reference_evaluation.json')],capture_output=True,text=True)
    (ROOT/'results/evaluator_output.txt').write_text(evaluation.stdout+evaluation.stderr)
    if evaluation.returncode:raise SystemExit(evaluation.returncode)
    command=[sys.executable,'-m','pytest','tests','-q','--junitxml=results/junit.xml']
    test=subprocess.run(command,capture_output=True,text=True,cwd=ROOT)
    (ROOT/'results/test_report.txt').write_text(test.stdout+test.stderr)
    print(test.stdout)
    if test.stderr:print(test.stderr,file=sys.stderr)
    if test.returncode:raise SystemExit(test.returncode)
    xml=ET.parse(ROOT/'results/junit.xml').getroot()
    suites=list(xml) if xml.tag=='testsuites' else [xml]
    summary={'tests':sum(int(s.get('tests',0)) for s in suites),'failures':sum(int(s.get('failures',0)) for s in suites),'errors':sum(int(s.get('errors',0)) for s in suites),'skipped':sum(int(s.get('skipped',0)) for s in suites),'reference_evaluator_score':json.loads((ROOT/'results/reference_evaluation.json').read_text())['score_out_of_100']}
    (ROOT/'results/validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()

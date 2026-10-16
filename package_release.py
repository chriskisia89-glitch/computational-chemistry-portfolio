"""Build clean candidate/instructor archives and a SHA-256 file manifest."""
from pathlib import Path
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parents[1]


def files():
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or any(x in p.parts for x in ['__pycache__','.pytest_cache','.venv']):continue
        if p.suffix=='.pyc' or p.name=='MANIFEST.json':continue
        yield p


def main():
    summary=json.loads((ROOT/'results/validation_summary.json').read_text())
    if summary['failures'] or summary['errors'] or summary['reference_evaluator_score']!=100:
        raise RuntimeError('release requires passing tests and reference evaluation')
    candidate=ROOT/'candidate_bundle.zip'
    with zipfile.ZipFile(candidate,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for source,name in [('task/TASK.md','TASK.md'),('task/isotherms.csv','isotherms.csv'),('task/submission_template.py','solution.py'),('requirements.txt','requirements.txt')]:
            z.write(ROOT/source,name)
        license_text=(ROOT/'LICENSE').read_text().split('This license covers')[0]
        z.writestr('LICENSE',license_text)
        z.writestr('README.md','# Candidate package\n\nRead TASK.md, complete solve(csv_path) in solution.py, then run:\n\n```bash\npython -m pip install -r requirements.txt\npython solution.py isotherms.csv answers.json\n```\n\nAll data in this candidate package are synthetic. No solution, ground truth, instructor validation dataset or source-derived plot data is included.\n')
    manifest={'release':'1.0.0','created_utc_date':'2026-10-16',
      'credit':'AI-assisted implementation; original experiments credited separately',
      'validation':summary,
      'files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files()]}
    (ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    final=ROOT.parent/'Chris_Kisia_Adsorption_Portfolio.zip'
    with zipfile.ZipFile(final,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in files():z.write(p,Path('adsorption_portfolio')/p.relative_to(ROOT))
        z.write(ROOT/'MANIFEST.json','adsorption_portfolio/MANIFEST.json')
    with zipfile.ZipFile(final) as z:
        if z.testzip() is not None:raise RuntimeError('archive CRC check failed')
        archive_names=set(z.namelist())
        assert not any('__pycache__' in n or '.pytest_cache' in n for n in archive_names)
        for record in manifest['files']:
            payload=z.read('adsorption_portfolio/'+record['path'])
            assert hashlib.sha256(payload).hexdigest()==record['sha256']
    with zipfile.ZipFile(candidate) as z:
        assert set(z.namelist())=={'TASK.md','isotherms.csv','solution.py','requirements.txt','LICENSE','README.md'}
    print(json.dumps({'archive':str(final),'bytes':final.stat().st_size,'files':len(archive_names),'validation':summary},indent=2))

if __name__=='__main__':main()

"""Reference implementation of the candidate-facing solve(csv_path) API."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from virial_portfolio.model import benchmark_answers


def solve(csv_path):
    return benchmark_answers(csv_path)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser();p.add_argument('csv_path');p.add_argument('output_json')
    a=p.parse_args();Path(a.output_json).write_text(json.dumps(solve(a.csv_path),indent=2,allow_nan=False)+'\n')

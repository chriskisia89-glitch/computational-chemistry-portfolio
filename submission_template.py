"""Complete solve(csv_path); input and output contracts are in TASK.md."""
from pathlib import Path


def solve(csv_path):
    raise NotImplementedError('Implement the specified two-temperature virial calculation')


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser();p.add_argument('csv_path');p.add_argument('output_json')
    a=p.parse_args();Path(a.output_json).write_text(json.dumps(solve(a.csv_path),indent=2,allow_nan=False)+'\n')

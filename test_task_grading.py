from copy import deepcopy
import json
from pathlib import Path
import pytest
from virial_portfolio.model import benchmark_answers
from grade_submission import validate_answer
from evaluate_solution import evaluate

ROOT=Path(__file__).resolve().parents[1]


def expected():return json.loads((ROOT/'task/expected_answers.json').read_text())


def test_reference_matches_frozen_gold():
    assert validate_answer(benchmark_answers(ROOT/'task/isotherms.csv'),expected())==[]


def test_grader_rejects_wrong_heat_units():
    a=expected();a['q0_kj_mol']*=1000
    assert validate_answer(a,expected())


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),True,'21'])
def test_grader_rejects_nonfinite_or_nonnumeric(bad):
    a=expected();a['q0_kj_mol']=bad
    assert validate_answer(a,expected())


def test_grader_rejects_coefficient_order_change():
    a=expected();a['coefficients']['a'].reverse()
    assert validate_answer(a,expected())


def test_grader_rejects_missing_keys_and_wrong_shapes():
    a=expected();a.pop('rank');assert validate_answer(a,expected())
    a=expected();a['coefficients']['b'].append(0);assert validate_answer(a,expected())


def test_grader_accepts_small_rounding_within_tolerance():
    a=expected();a['q0_kj_mol']+=1e-8
    assert validate_answer(a,expected())==[]


def test_second_dataset_matches_independent_qr_reference():
    gold=json.loads((ROOT/'tests/fixtures/validation_expected.json').read_text())
    assert validate_answer(benchmark_answers(ROOT/'tests/fixtures/validation_isotherms.csv'),gold)==[]


def test_reference_passes_end_to_end_evaluator():
    result=evaluate(ROOT/'task/reference_solution.py')
    assert result['passed_all_cases'] and result['score_out_of_100']==100


def test_dataset_lookup_is_not_a_complete_solution(tmp_path):
    solution=tmp_path/'lookup.py'
    solution.write_text('import json\nfrom pathlib import Path\ndef solve(path):\n    return json.loads(Path('+repr(str(ROOT/'task/expected_answers.json'))+').read_text())\n')
    result=evaluate(solution)
    assert not result['passed_all_cases'] and result['score_out_of_100']<100


def test_evaluator_supports_candidate_dataclasses(tmp_path):
    solution=tmp_path/'dataclass_solution.py'
    solution.write_text('from dataclasses import dataclass\nfrom virial_portfolio.model import benchmark_answers\n@dataclass\nclass Options:\n    enabled: bool = True\ndef solve(path):\n    return benchmark_answers(path) if Options().enabled else {}\n')
    assert evaluate(solution)['score_out_of_100']==100


def test_evaluator_awards_correct_component_partial_credit(tmp_path):
    solution=tmp_path/'wrong_heat_units.py'
    solution.write_text('from virial_portfolio.model import benchmark_answers\ndef solve(path):\n    result=benchmark_answers(path)\n    result["q0_kj_mol"]*=1000\n    return result\n')
    result=evaluate(solution)
    assert result['score_out_of_100']==80
    assert not result['passed_all_cases']

"""Production benchmark builders and runner (M0): must reproduce the prototype benchmarks and enforce SPEC 7.4 floors."""
import json, os, pickle, subprocess, sys
import pytest
from rdkit import Chem

from open_cgenff.benchmarks import b0, b1, b3, floors, runner


@pytest.fixture(scope='session')
def b1_cases(toppar):
    return b1.build(os.path.join(toppar, 'top_all36_cgenff.rtf'))


@pytest.fixture(scope='session')
def report(toppar):
    return runner.run(toppar=toppar)


def test_b1_matches_prototype_build(b1_cases, b1):
    cases, stats = b1_cases
    assert stats['molecules'] == 1130 and stats['atoms'] == 27512 and stats['repaired_radicals'] == 94
    proto = {d['name']: d for d in b1}
    assert set(proto) == {c.name for c in cases}
    for c in cases:
        d = proto[c.name]
        assert c.types == d['types'] and c.split == d['split'], c.name
        assert Chem.MolToSmiles(c.mol) == Chem.MolToSmiles(Chem.Mol(d['mol'])), c.name


def test_b1_manifest_is_reproducible(b1_cases, root):
    manifest = json.load(open(os.path.join(root, 'benchmarks', 'b1_rtf_zinc.split.json')))
    assert b1.shuffle_split(b1_cases[0]) == {'train': manifest['train'], 'test': manifest['test']}


def test_b3_matches_prototype_build(root, tmp_path):
    out = tmp_path / 'b3.pkl'
    subprocess.run([sys.executable, os.path.join(root, 'prototype', 'charmmgui_benchmark.py'),
                    os.path.join(root, 'benchmarks', 'charmm_gui_cpa', 'raw'),
                    os.path.join(root, 'prototype', 'typer_pilot.py'), '--out', str(out)], check=True, capture_output=True)
    proto = {d['file'].removeprefix('cgen_'): d for d in pickle.load(open(out, 'rb'))}
    cases, excluded = b3.build()
    assert len(cases) == 34 and len(excluded) == 5
    for c in cases:
        d = proto[c.meta['file'].removeprefix('cgen_')]
        assert c.types == d['types'] and c.meta['charges'] == d['charges']
        assert Chem.MolToSmiles(c.mol) == Chem.MolToSmiles(Chem.Mol(d['mol']))


def test_runner_meets_floors(report):
    assert report['floor_violations'] == []
    assert report['B0']['exact'] == 3
    assert report['B1']['typing']['train']['molecules_correct'] >= 874
    assert report['B3']['typing']['dev']['molecule_accuracy'] == 1.0


def test_runner_hides_sealed_split_details(report):
    test = report['B1']['typing']['test']
    assert 'confusions' not in test and 'failures' not in test
    assert 'confusions' in report['B1']['typing']['train']
    assert 'Top confusions, B1-test' not in runner.to_markdown(report)


def test_floor_check_detects_regression():
    bad = {'B1': {'typing': {'train': {'molecules_correct': 873, 'molecules': 904},
                             'test': {'molecules_correct': 212, 'molecules': 226}}},
           'B0': {'exact': 2, 'molecules': 3}}
    assert len(floors.check(bad)) == 2


def test_runner_cli_writes_reports(toppar, tmp_path):
    rc = runner.main(['--sets', 'B0,B3', '--toppar', toppar, '--out', str(tmp_path), '--check-floors'])
    assert rc == 0
    data = json.load(open(tmp_path / 'report.json'))
    assert 'not the official CGenFF program' in data['tool']
    assert (tmp_path / 'report.md').read_text().startswith('# Benchmark report')

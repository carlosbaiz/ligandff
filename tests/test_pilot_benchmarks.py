"""Regression floor for the pilot typer. Production code must meet or beat these numbers (docs/SPEC.md section 7)."""
import os, pickle, subprocess, sys, importlib
import pytest
from rdkit import Chem

ALT = {'CG2DC1': 'CG2DC?', 'CG2DC2': 'CG2DC?', 'CG2D1O': 'CG2D?O', 'CG2D2O': 'CG2D?O',
       'CG25C1': 'CG25C?', 'CG25C2': 'CG25C?', 'CG251O': 'CG25?O', 'CG252O': 'CG25?O'}


def molecule_accuracy(data, type_molecule, split=None):
    rows = [d for d in data if split is None or d['split'] == split]
    ok = sum(all(ALT.get(t, t) == ALT.get(p, p) for t, p in zip(d['types'], type_molecule(Chem.Mol(d['mol'])))) for d in rows)
    return ok / len(rows)


@pytest.fixture(scope='session')
def b1(toppar, root, tmp_path_factory):
    out = tmp_path_factory.mktemp('b1') / 'b1.pkl'
    subprocess.run([sys.executable, os.path.join(root, 'prototype', 'build_rtf_benchmark.py'),
                    os.path.join(toppar, 'top_all36_cgenff.rtf'), str(out)], check=True, capture_output=True)
    return pickle.load(open(out, 'rb'))


def test_b1_split_is_stable(b1, root):
    import json
    ref = json.load(open(os.path.join(root, 'benchmarks', 'b1_rtf_zinc.split.json')))
    assert sorted(d['name'] for d in b1 if d['split'] == 'test') == ref['test']


def test_pilot_typer_b1(b1):
    typer = importlib.import_module('typer_pilot')
    assert molecule_accuracy(b1, typer.type_molecule, 'train') >= 0.96
    assert molecule_accuracy(b1, typer.type_molecule, 'test') >= 0.93


def test_pilot_typer_charmmgui(root, tmp_path):
    out = tmp_path / 'b3.pkl'
    subprocess.run([sys.executable, os.path.join(root, 'prototype', 'charmmgui_benchmark.py'),
                    os.path.join(root, 'benchmarks', 'charmm_gui_cpa', 'raw'),
                    os.path.join(root, 'prototype', 'typer_pilot.py'), '--out', str(out)], check=True, capture_output=True)
    data = pickle.load(open(out, 'rb'))
    typer = importlib.import_module('typer_pilot')
    assert len(data) == 34
    assert molecule_accuracy(data, typer.type_molecule) == 1.0

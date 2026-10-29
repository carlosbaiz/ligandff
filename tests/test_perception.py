"""Perception (M2, SPEC 5.1-5.4): mol2 input, rings, resonance (P-I worked examples), aromaticity, special cases."""
import glob, os
import pytest
from rdkit import Chem

from open_cgenff.io.mol2 import read_mol2
from open_cgenff.molecule import AROMATIC as AR, InputMolecule
from open_cgenff.perception import COEFFICIENTS_PRE_097, find_rings, perceive, resolve

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECIAL = os.path.join(ROOT, 'data', 'pi_si', 'specialcases')


def build(heavy, edges, n_h, order=AR):
    """Heavy atoms 0..k-1 joined by `edges` (bond order `order` unless given as (a, b, o)), then n_h[i] hydrogens."""
    el = list(heavy)
    bonds = [(e[0], e[1], e[2] if len(e) == 3 else order) for e in edges]
    for i, k in enumerate(n_h):
        for _ in range(k):
            el.append('H'); bonds.append((i, len(el) - 1, 1))
    return InputMolecule('x', [f'A{i}' for i in range(len(el))], el, bonds)


def ring(n):
    return [(i, (i + 1) % n) for i in range(n)]


def smiles(s):
    return perceive(InputMolecule.from_rdkit(Chem.AddHs(Chem.MolFromSmiles(s))))


# ---- P-I worked examples ---------------------------------------------------------------------------------------

def fig4_molecule():
    edges = [(1, 2), (2, 3), (3, 4), (3, 10), (4, 5), (4, 9), (5, 6), (6, 7), (7, 8), (8, 9), (7, 11), (11, 12)]
    return build('CCCCCNCCCCCO', [(a - 1, b - 1) for a, b in edges], [2, 1, 0, 0, 1, 1, 0, 1, 1, 2, 1, 0])


def test_fig4_structures_found_in_published_order():
    """P-I Fig. 4: valid structures H (12), I (11), J (11) in that order; I (N6+, N6=C7, C11=O12) is kept."""
    p = perceive(fig4_molecule())
    assert p.resonance.structures == [12, 11, 11]
    assert p.resonance.score == 11
    assert p.charges[:12] == [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    assert (p.order(5, 6), p.order(10, 11), p.order(3, 4), p.order(7, 8), p.order(9, 2)) == (2, 2, 2, 2, 2)
    assert list(p.ring_class.values()) == ['arom']


@pytest.mark.parametrize('coefficients, score', [(None, 11), (COEFFICIENTS_PRE_097, 6)])
def test_pyridinium_two_structures(coefficients, score):
    """P-I case study 1: two structures, score 6 with the pre-0.9.7 coefficients (4 + 2) and 8 + 3 = 11 now."""
    m = build('NCCCCC', ring(6), [1] * 6)
    nb = [sorted(j for i2, j2, _ in m.bonds for i, j in ((i2, j2), (j2, i2)) if i == a) for a in range(m.n_atoms)]
    kw = {'coefficients': coefficients} if coefficients else {}
    r = resolve(m.n_atoms, m.elements, nb, m.bonds, None, find_rings(nb), keep_all=True, **kw)
    assert [s[0] for s in r.structures] == [score, score] and r.charges[0] == 1


def test_benzene_stops_at_score_zero():
    p = perceive(build('CCCCCC', ring(6), [1] * 6))
    assert p.resonance.structures == [0] and p.order(0, 1) == 2 and list(p.ring_class.values()) == ['arom']


def test_pyrrole_nitrogen_neutral():
    p = perceive(build('NCCCC', ring(5), [1] * 5))
    assert p.charges[0] == 0 and p.resonance.score == 0 and list(p.ring_class.values()) == ['arom']


def test_charges_outside_the_subgraph_are_not_scored():
    """A fixed dianion elsewhere must not make an N+ (non-aromatic) structure win (FINDINGS 6)."""
    p = smiles('c1ccncc1CP(=O)([O-])[O-]')
    assert all(c == 'arom' for c in p.ring_class.values())
    assert sum(p.charges) == -2 and p.resonance.score == 0


def test_carboxylate_with_aromatic_bonds():
    p = perceive(build('CCOO', [(0, 1, 1), (1, 2), (1, 3)], [3, 0, 0, 0]))
    assert sorted(p.charges[2:4]) == [-1, 0] and sorted((p.order(1, 2), p.order(1, 3))) == [1, 2]


def test_carbanion_fallback_warns():
    p = perceive(build('CCCCC', ring(5), [1] * 5))          # cyclopentadienide given with aromatic bonds
    assert any('carbanion' in w for w in p.warnings) and sum(p.charges) == -1
    assert list(p.ring_class.values()) == ['arom']


# ---- rings -----------------------------------------------------------------------------------------------------

def test_rings_naphthalene():
    p = smiles('c1ccc2ccccc2c1')
    assert sorted(map(len, p.rings)) == [6, 6]                  # the 10-ring is too large
    fusion = [i for i in range(10) if p.bridgehead[i]]
    assert len(fusion) == 2 and all(len(p.atom_rings[i]) == 2 for i in fusion)
    assert all(c == 'arom' for c in p.ring_class.values())


def test_rings_norbornane():
    p = smiles('C1CC2CCC1C2')
    assert sorted(map(len, p.rings)) == [5, 5, 6]
    heads = [i for i in range(7) if p.bridgehead[i]]
    assert len(heads) == 2 and all(sorted(map(len, p.atom_rings[i])) == [5, 5, 6] for i in heads)
    assert all(len(p.atom_rings[i]) == 1 and len(p.atom_rings[i][0]) == 5 for i in range(7) if i not in heads)


def test_rings_bicyclooctane_spiro_cubane_cyclooctane():
    assert sorted(map(len, smiles('C1CC2CCC1CC2').rings)) == [6, 6, 6]
    spiro = smiles('C1CCC2(C1)CCCC2')
    assert sorted(map(len, spiro.rings)) == [5, 5]
    cubane = smiles('C12C3C4C1C5C2C3C45')
    assert sum(len(r) == 4 for r in cubane.rings) == 6
    assert all(sorted(map(len, cubane.atom_rings[i])) == [4, 4, 4] for i in range(8))
    assert smiles('C1CCCCCCC1').rings == []


# ---- ring classes ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize('smi, classes', [
    ('c1ccncc1', ['arom']), ('c1ccoc1', ['arom']), ('O=c1cccc[nH]1', ['arom']), ('C1=CCCCC1', ['mixed']),
    ('C1CCCCC1', ['sp3']), ('C1=CCC=C1', ['mixed']), ('c1ccc2[nH]ccc2c1', ['arom', 'arom']),
    ('C1=CC=CC1=O', ['sp2']), ('C1CC1', ['sp3']), ('c1ccc2occc2c1', ['arom', 'arom']),
])
def test_ring_classes(smi, classes):
    p = smiles(smi)
    assert sorted(p.ring_class[r] for r in p.rings) == sorted(classes)


# ---- mol2 input and P-I special cases --------------------------------------------------------------------------

def test_mol2_reader(tmp_path):
    f = tmp_path / 'pyr.mol2'
    atoms = ''.join(f'{i + 1} {e}{i + 1} 0 0 0 {t}\n' for i, (e, t) in enumerate(
        [('N', 'N.ar'), ('C', 'C.ar'), ('C', 'C.ar'), ('C', 'C.ar'), ('C', 'C.ar'), ('C', 'C.ar'), ('H', 'H')]))
    bonds = ''.join(f'{k + 1} {a + 1} {b + 1} ar\n' for k, (a, b) in enumerate(ring(6))) + '7 1 7 1\n'
    f.write_text(f'@<TRIPOS>MOLECULE\npyridinium\n7 7\nSMALL\n\n@<TRIPOS>ATOM\n{atoms}@<TRIPOS>BOND\n{bonds}')
    m = read_mol2(f)[0]
    assert m.name == 'pyridinium' and m.elements[:2] == ['N', 'C'] and m.formal_charges is None
    assert sum(o == AR for _, _, o in m.bonds) == 6
    p = perceive(m)
    assert p.charges[0] == 1


@pytest.fixture(scope='session')
def special_cases():
    files = sorted(glob.glob(os.path.join(SPECIAL, '*.mol2')))
    if len(files) != 126:            # copyrighted SI: imported locally with scripts/import_pi_si.sh, absent in CI
        pytest.skip('P-I SI special cases not imported (scripts/import_pi_si.sh)')
    return {os.path.basename(f)[:-5]: read_mol2(f)[0] for f in files}


def test_all_special_cases_run(special_cases):
    warned = {}
    for name, m in special_cases.items():
        p = perceive(m)
        assert len(p.charges) == m.n_atoms
        if p.warnings:
            warned[name] = p.warnings
    assert set(warned) == {'fuller_all_ar', 'isoindenate'}, warned


def test_special_cases_quinolinium_and_spiro(special_cases):
    p = perceive(special_cases['quinolinium'])
    assert sum(p.charges) == 1 and all(c == 'arom' for c in p.ring_class.values())
    s = perceive(special_cases['spiro1'])
    assert s.rings and all(c == 'sp3' for c in s.ring_class.values())


# ---- benchmarks: B2 and the pilot typer on production perception -------------------------------------------------

def test_b2_build(toppar):
    from open_cgenff.benchmarks import b1, b2
    cases, stats = b2.build(os.path.join(toppar, 'top_all36_cgenff.rtf'))
    assert stats['molecules'] == 910 and {c.split for c in cases} == {'train', 'test'}
    b1_types = {t for c in b1.build(os.path.join(toppar, 'top_all36_cgenff.rtf'))[0] for t in c.types}
    assert len({t for c in cases for t in c.types} - b1_types) == 32


def test_pilot_typer_on_production_perception(toppar):
    """SPEC 5.3: pilot results with production perception (B1-train -1 / B1-test +1 molecule vs the pilot's own
    perception, all explained in FINDINGS 6); aggregate counts only for test splits."""
    from open_cgenff.benchmarks import runner
    report = runner.run(['B1', 'B2', 'B3'], typer=os.path.join(ROOT, 'benchmarks', 'plugins',
                                                               'pilot_on_production_perception.py'), toppar=toppar)
    got = {(s, sp): m['molecules_correct'] for s in ('B1', 'B2', 'B3') for sp, m in report[s]['typing'].items()}
    assert got[('B1', 'train')] >= 873 and got[('B1', 'test')] >= 213
    assert got[('B2', 'train')] >= 627 and got[('B2', 'test')] >= 158 and got[('B3', 'dev')] == 34

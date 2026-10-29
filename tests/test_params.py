"""Parameter database, penalty rules and the analogy engine (SPEC 5.8-5.10, M4)."""
import json
import os

import pytest

from open_cgenff.charges import IncrementTable, charge_molecule
from open_cgenff.molecule import TypedMolecule
from open_cgenff.params import (AnalogyEngine, BondGroups, ParameterDatabase, PenaltyModel, PenaltyRuleError,
                                PenaltyRules, SubstitutionMatrix, assign_parameters, charge_analogs, default_engine,
                                parse_penalty_rules)
from open_cgenff.params.penalties import BondGroup

FIXTURES = os.path.join(os.path.dirname(__file__), 'fixtures')

# ---- synthetic database and model -----------------------------------------------------------------------------------

TOY_PRM = """\
ATOMS
MASS  -1  CA   12.011
MASS  -1  CB   12.011
MASS  -1  CC   12.011
MASS  -1  NA   14.007
MASS  -1  HA    1.008
BONDS
CA   CA   300.0  1.40
CA   NA   200.0  1.30
ANGLES
CA   CA   CA    40.0  120.0
NA   CA   NA    50.0  110.0   10.0  2.20
DIHEDRALS
CA   CA   CA   CA     1.0  2  180.0
HA   CA   CA   NA     0.5  3    0.0
HA   CA   CA   NA     0.2  1    0.0
IMPROPERS
CA   NA   NA   HA    40.0  0    0.0
END
"""


@pytest.fixture
def toy_db(tmp_path):
    path = tmp_path / 'toy.prm'
    path.write_text(TOY_PRM)
    return ParameterDatabase.from_prm(path)


def toy_model(bonded=None, nonbonded=None, groups=()):
    """Matrices with given elements; every other off-diagonal element costs 100."""
    fallback = {r: 100.0 for r in (('diff',), ('elem', False), ('elem', True), ('fam', 0), ('fam', 1), ('fam', 2))}
    return PenaltyModel(SubstitutionMatrix('bonded', bonded or {}, fallback),
                        SubstitutionMatrix('nonbonded', nonbonded or {}, fallback),
                        BondGroups(groups))


def test_toy_database_reads_all_terms(toy_db):
    assert len(toy_db.params['bond']) == 2 and len(toy_db.params['improper']) == 1
    ang = toy_db.params['angle'][('NA', 'CA', 'NA')]
    assert ang.values == (50.0, 110.0, 10.0, 2.2)                 # Urey-Bradley term
    dih = toy_db.params['dihedral'][('HA', 'CA', 'CA', 'NA')]
    assert dih.values == ((0.5, 3, 0.0), (0.2, 1, 0.0))            # multi-term dihedral
    assert toy_db.lookup('bond', ('NA', 'CA'))[1] == -1


def test_inner_positions_weigh_ten():
    eng = AnalogyEngine(None, toy_model(bonded={('CB', 'CA'): 1.0}, nonbonded={('CB', 'CA'): 1.0}))
    assert eng.tps('bond', ('CB', 'CA'), ('CA', 'CA')) == 10.0
    assert eng.tps('angle', ('CB', 'CA', 'CA'), ('CA', 'CA', 'CA')) == 1.0     # end: nonbonded, x1
    assert eng.tps('angle', ('CA', 'CB', 'CA'), ('CA', 'CA', 'CA')) == 10.0    # centre: bonded, x10
    assert eng.tps('dihedral', ('CA', 'CB', 'CA', 'CA'), ('CA',) * 4) == 10.0
    assert eng.tps('dihedral', ('CA', 'CA', 'CA', 'CB'), ('CA',) * 4) == 1.0
    assert eng.tps('improper', ('CB', 'CA', 'CA', 'CA'), ('CA',) * 4) == 10.0  # centre x10
    assert eng.tps('improper', ('CA', 'CA', 'CB', 'CA'), ('CA',) * 4) == 1.0   # outer: bonded matrix, x1


def test_matrices_are_used_by_position_and_are_asymmetric():
    eng = AnalogyEngine(None, toy_model(bonded={('CB', 'CA'): 1.0}, nonbonded={('CB', 'CA'): 3.0, ('CA', 'CB'): 7.0}))
    assert eng.tps('angle', ('CB', 'CA', 'CA'), ('CA', 'CA', 'CA')) == 3.0
    assert eng.tps('angle', ('CA', 'CA', 'CA'), ('CB', 'CA', 'CA')) == 7.0
    # charge-increment analogs: nonbonded matrix at every position, same weights (SPEC 5.11)
    assert eng.tps('angle', ('CA', 'CB', 'CA'), ('CA', 'CA', 'CA'), mode='charge') == 30.0


def test_bond_group_penalties():
    groups = [BondGroup(40, frozenset({'CA', 'CB'})), BondGroup(40, frozenset({'CC'})),   # equivalent pair
              BondGroup(20, frozenset({'CA', 'NA'}))]
    sub = {('CB', 'CC'): 0.0, ('NA', 'CC'): 0.0, ('CC', 'NA'): 0.0}
    eng = AnalogyEngine(None, toy_model(bonded=sub, nonbonded=sub, groups=groups))
    # bond parameter: virtual bond weighs 10; CA-CB in group 1, CA-CC in no group
    assert eng.tps('bond', ('CA', 'CB'), ('CA', 'CC')) == 400.0
    # first two groups are equivalent: CB-CB (group 1) vs CC-CC (group 2) is no difference
    assert eng.tps('bond', ('CB', 'CB'), ('CC', 'CC')) == 0.0
    # dihedral: outer bonds x1, central bond x10; CA-NA (group 3) vs CA-CC (none) at an outer bond
    assert eng.tps('dihedral', ('NA', 'CA', 'CA', 'CA'), ('CC', 'CA', 'CA', 'CA')) == 20.0
    assert eng.tps('dihedral', ('CA', 'NA', 'CA', 'CA'), ('CA', 'CC', 'CA', 'CA')) == 20.0 + 200.0
    # impropers have no virtual bonds
    assert eng.tps('improper', ('NA', 'CA', 'CA', 'CA'), ('CC', 'CA', 'CA', 'CA')) == 0.0


def test_search_reversal_ties_and_cache(toy_db):
    model = toy_model(bonded={('CB', 'CA'): 0.5, ('CC', 'CA'): 2.0, ('CB', 'NA'): 0.5},
                      nonbonded={('HA', 'NA'): 1.0, ('NA', 'HA'): 1.0, ('CB', 'CA'): 0.2})
    eng = AnalogyEngine(toy_db, model)
    # NA-CB: as given would need NA->CA (100); reversed CB-NA vs CA-NA costs 10 x 0.5
    res = eng.find('bond', ('NA', 'CB'))
    assert res.source == ('CA', 'NA') and res.penalty == 5.0 and res.aligned == ('NA', 'CA')
    # tie: CB-CB vs CA-CA (10 x 0.5 x 2) or CA-NA (100s): CA-CA; equal candidates -> first in collated order
    res = eng.find('angle', ('CA', 'CA', 'CB'))
    assert res.source == ('CA', 'CA', 'CA') and res.penalty == 0.2
    assert eng.find('angle', ('CA', 'CA', 'CB')) is res                  # cached
    # outer substitution at a dihedral end (nonbonded matrix)
    res = eng.find('dihedral', ('NA', 'CA', 'CA', 'NA'))
    assert res.source == ('HA', 'CA', 'CA', 'NA') and res.penalty == 1.0
    # improper: outer-atom orders are tried; aligned follows the missing order
    res = eng.find('improper', ('CA', 'HA', 'NA', 'NA'))
    assert res.penalty == 0.0 and res.aligned == ('CA', 'HA', 'NA', 'NA')


def test_collated_order_breaks_ties(tmp_path):
    prm = TOY_PRM.replace('CA   NA   200.0  1.30', 'CA   NA   200.0  1.30\nCA   CC   250.0  1.35')
    (tmp_path / 'p.prm').write_text(prm)
    db = ParameterDatabase.from_prm(tmp_path / 'p.prm')
    eng = AnalogyEngine(db, toy_model(bonded={('CB', 'CC'): 1.0, ('CB', 'NA'): 1.0}))
    assert [p.types for p in db.collated('bond')] == [('CA', 'CA'), ('CA', 'CC'), ('CA', 'NA')]
    assert eng.find('bond', ('CA', 'CB')).source == ('CA', 'CC')


# ---- penalty rule file ------------------------------------------------------------------------------------------------

def test_published_ng3_hierarchy():
    rules = PenaltyRules.from_file()
    cost = rules.hierarchy_cost
    assert cost('bonded', 'NG3P3', 'NG321') == 10.0         # P-II worked example: 8 + 2 + 0
    assert cost('bonded', 'NG3P3', 'NG311') == 10.5         # 8 + 2 + 0.5
    assert cost('bonded', 'NG311', 'NG301') == 0.5 and cost('bonded', 'NG301', 'NG311') == 0.5
    assert cost('bonded', 'NG3C51', 'NG311') == 2.0
    assert cost('nonbonded', 'NG3P3', 'NG321') is None      # Table 1 is the bonded matrix
    assert cost('bonded', 'NG321', 'CG321') is None
    assert len(rules.bond_groups) == 10 and 'CG2RC7' in rules.bond_groups[5].types   # continuation line


@pytest.mark.parametrize('text', [
    'cat A\ntyp X : pri 0 alt Y 1\nend\n',                  # alt to an entry not in the category
    'cat A\nsub B : pri 0\nend\n',                          # sub to an undefined category
    'typ X : pri 0\n',                                      # entry outside a category
    'cat A\ntyp X : pri zero\nend\n',
    'matrix other\n',
    'frob\n',
])
def test_penalty_rule_errors(text):
    with pytest.raises(PenaltyRuleError):
        parse_penalty_rules(text)


# ---- real data ---------------------------------------------------------------------------------------------------------

@pytest.fixture(scope='module')
def engine(toppar):
    return default_engine(toppar)


def test_parameter_file_and_collated_order(engine):
    db = engine.db
    assert {k: len(v) for k, v in db.params.items()} == {'bond': 1305, 'angle': 6776, 'dihedral': 19629,
                                                         'improper': 489}
    assert sum(p.annotation is not None for t in db.params.values() for p in t.values()) == 7234
    for kind in ('bond', 'angle', 'dihedral'):
        assert all(db.canonical(kind, p.types) == p.types for p in db.params[kind].values())
    for kind in ('bond', 'angle'):                          # the file is written in collated order
        assert db.collated(kind) == list(db.params[kind].values())
    assert db.elements['CLGR1'][1] == 'Cl' and db.priority('CG2O5') < db.priority('CG311') < db.priority('HGA1')


def test_provisional_matrix_sources(engine):
    m = engine.model
    assert m.bonded.source('NG3P3', 'NG321') == 'rules' and m.bonded.cost('NG3P3', 'NG321') == 10.0
    assert m.nonbonded.source('CG311', 'CG321') == 'observed'
    assert m.bonded.source('CG2R61', 'SG311') == 'fallback'
    assert m.bonded.cost('CG2R61', 'CG2R61') == 0.0


def load_fixture(name):
    fx = {m['name']: m for m in json.load(open(os.path.join(FIXTURES, 'si_v5_reference.json')))}[name]
    names = [a['name'] for a in fx['atoms']]
    idx = {n: i for i, n in enumerate(names)}
    mol = TypedMolecule(names, [a['type'] for a in fx['atoms']], [(idx[a], idx[b]) for a, b in fx['bonds']],
                        [fx['formal_charges'].get(n, 0) for n in names])
    return fx, mol


@pytest.mark.parametrize('name', ['M064', 'M166', 'M223'])
def test_si_fixture_analogs_penalties_and_charges(engine, toppar, name):
    """M4 acceptance: the program 5.0 analog of M064 (SI Table S10) is reproduced, and with analog increments the
    charges and charge penalties match the official output."""
    fx, mol = load_fixture(name)
    asg = assign_parameters(mol, engine)
    assert not asg.missing
    got = [{'missing': list(k[1]), 'source': list(a.source), 'penalty': a.penalty} for k, a in asg.analogs().items()]
    assert got == fx['analog_parameters']
    assert asg.param_penalty == fx['param_penalty']
    table = IncrementTable.from_toppar(toppar)
    r = charge_molecule(mol, table, charge_analogs(mol, engine, table))
    assert max(abs(q - a['charge']) for q, a in zip(r.charges, fx['atoms'])) < 5e-4
    assert max(abs(p - a['charge_penalty']) for p, a in zip(r.penalties, fx['atoms'])) < 5e-4
    assert round(r.charge_penalty, 3) == fx['charge_penalty']


def test_typed_molecule_gets_improper_parameters(engine):
    from rdkit import Chem
    from open_cgenff.molecule import InputMolecule
    from open_cgenff.typing import type_molecule
    tm = type_molecule(InputMolecule.from_rdkit(Chem.AddHs(Chem.MolFromSmiles('CC(=O)NC'))))
    asg = assign_parameters(tm, engine)
    imps = [p for p in asg.parameters if p.kind == 'improper']
    assert len(imps) == 1 and imps[0].penalty == 0.0 and imps[0].types[0] == 'CG2O1'
    assert tm.types[imps[0].atoms[0]] == 'CG2O1'
    assert asg.param_penalty == 0.0 and not asg.missing


def brute_force(engine, kind, types, mode='bonded', accept=None):
    """Reference search: every candidate in collated order, every orientation; first lowest TPS."""
    best = None
    for prm in engine.db.collated(kind):
        if accept is not None and not accept(prm):
            continue
        for perm in engine.orientations(kind, len(types)):
            score = engine.tps(kind, tuple(types[i] for i in perm), prm.types, mode)
            if best is None or score < best[0] - 1e-9:
                best = (score, prm)
    return best


def test_pruned_search_equals_brute_force(engine, toppar):
    """Leave-one-out on annotated parameters of each kind: the indexed search returns what a full scan returns."""
    from open_cgenff.params.analogy import _increment_filter
    table = IncrementTable.from_toppar(toppar)
    picks = {'bond': 2, 'angle': 2, 'dihedral': 2, 'improper': 2}
    for kind, n in picks.items():
        annotated = [p for p in engine.db.params[kind].values() if p.annotation][::97][:n]
        for prm in annotated:
            db = engine.db.without(kind, prm.types)
            eng = AnalogyEngine(db, engine.model)
            for mode, accept in (('bonded', None), ('charge', _increment_filter(table))):
                if kind == 'improper' and mode == 'charge':
                    continue
                res = eng.find(kind, prm.types, mode, accept)
                score, src = brute_force(eng, kind, prm.types, mode, accept)
                assert (res.source, res.penalty) == (src.types, round(score, 6)), (kind, prm.types, mode)


# ---- M5: fitted model and B4 --------------------------------------------------------------------------------------

def test_b4_split_is_reproducible(engine, root):
    from open_cgenff.params import annotation_split
    manifest = json.load(open(os.path.join(root, 'benchmarks', 'b4_annotations.split.json')))
    assert annotation_split(engine.db) == manifest
    assert len(manifest['train']) == 5787 and len(manifest['test']) == 1447


def test_fitted_model_conventions(engine):
    from open_cgenff.params.penalties import ALTNUM_TWINS, FALLBACK_SCALE
    m = engine.model
    for a, b in ALTNUM_TWINS:
        assert m.bonded.cost(a, b) == 0.0 and m.nonbonded.cost(b, a) == 0.0
    # unobserved substitutions are expensive (relation median x FALLBACK_SCALE)
    assert FALLBACK_SCALE == 20.0 and min(m.nonbonded.fallback.values()) >= 20.0
    assert m.bonded.cost('NG3P3', 'NG321') == 10.0                       # P-II Table 1 still wins


def test_b4_pool_excludes_the_batch_era(engine):
    from open_cgenff.benchmarks import b4
    cases = {c.name: c for c in b4.build(engine.db)}
    batch = next(c for c in cases.values() if c.batch and c.split == 'train')
    res = b4.leave_one_out(engine, batch)
    assert res is not None and not res.parameter.comment.startswith('aa,') and res.source != batch.types

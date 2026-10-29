"""Production charge engine (M1, SPEC 5.11): B0 exact, training check = FINDINGS 2.2, identical to the prototype."""
import math, os, random
import pytest

import charges as proto                                  # prototype/charges.py (reference baseline)
from open_cgenff.benchmarks import b0
from open_cgenff.charges import (Analog, IncrementTable, assign_charges, atom_penalty, charge_molecule,
                                 parameter_instances, training_check)
from open_cgenff.io.prm import read_increments
from open_cgenff.io.rtf import read_rtf
from open_cgenff.molecule import TypedMolecule


@pytest.fixture(scope='session')
def table(toppar):
    return IncrementTable.from_toppar(toppar)


def test_b0_exact_via_typed_molecule(table):
    for mol in b0.load():
        types, bonds, q0, analogs = b0.engine_inputs(mol)
        tm = TypedMolecule([a['name'] for a in mol['atoms']], types, bonds, q0)
        r = charge_molecule(tm, table, {k: Analog(*v) for k, v in analogs.items()})
        assert not r.missing
        assert [round(x, 3) for x in r.charges] == [a['charge'] for a in mol['atoms']], mol['name']
        assert [round(x, 3) for x in r.penalties] == [a['charge_penalty'] for a in mol['atoms']], mol['name']
        assert round(r.charge_penalty, 3) == mol['charge_penalty']
        assert set(r.analogs_used) <= set(analogs) | {k[::-1] for k in analogs}
        assert math.isclose(sum(r.charges), tm.net_charge, abs_tol=1e-9)


def test_training_check_reproduces_findings(toppar, table):
    r = training_check(toppar, table)
    assert r.residues == 2255
    assert round(r.median_rmse, 4) == 0.0098
    assert round(r.fraction_within_001, 3) == 0.760 and round(r.fraction_within_005, 3) == 0.940
    # FINDINGS 2.2 / SPEC 5.7: carboxylate C -1, ammonium N +1, phosphate P -1/-2
    for key in [('CG2O3', -1), ('NG3P2', 1), ('NG3P3', 1), ('PG1', -1), ('PG2', -2), ('CG2N2', 1)]:
        assert r.formal_charges.get(key, 0) > 0, key


def _training_residues(toppar):
    for res in read_rtf(os.path.join(toppar, 'top_all36_cgenff.rtf')):
        if res.kind != 'RESI' or res.has_external_bond or res.has_lonepair or len(res.atoms) < 2:
            continue
        idx = {a.name: i for i, a in enumerate(res.atoms)}
        bonds = sorted({tuple(sorted((idx[a], idx[b]))) for a, b in res.bonds if a in idx and b in idx})
        yield res, [a.type for a in res.atoms], bonds


def test_identical_to_prototype_on_training_residues(toppar, table):
    ptable = proto.load_increments(os.path.join(toppar, 'par_all36_cgenff_increments.prm'))
    n = 0
    for res, types, bonds in _training_residues(toppar):
        q0 = [round(a.charge) for a in res.atoms]
        r = assign_charges(types, bonds, q0, table)
        pq, ppen, pmissing = proto.assign(types, bonds, q0, ptable)
        assert r.charges == pq and r.penalties == ppen and r.missing == pmissing, res.name
        n += 1
    assert n > 2000


def test_identical_to_prototype_with_analogs(toppar, table):
    """Drop random parameters from residues and supply random analogs (both orientations, penalties above and below
    the dihedral cut-off): charges, penalties and missing lists must equal the prototype's."""
    rng = random.Random(7)
    ptable = proto.load_increments(os.path.join(toppar, 'par_all36_cgenff_increments.prm'))
    raw = read_increments(os.path.join(toppar, 'par_all36_cgenff_increments.prm'))
    residues = [x for x in _training_residues(toppar) if len(x[1]) >= 6]
    pool = {n: [k for k, v in raw[n].items() if v is not None] for n in (2, 3, 4)}
    for res, types, bonds in rng.sample(residues, 150):
        keys = {tuple(types[a] for a in inst) for inst in parameter_instances(len(types), bonds)}
        dropped = rng.sample(sorted(keys), min(4, len(keys)))
        sub_raw = {n: {k: v for k, v in raw[n].items() if k not in dropped and k[::-1] not in dropped} for n in raw}
        sub_proto = {n: {k: (list(v) if v else v) for k, v in d.items()} for n, d in sub_raw.items()}
        analogs = {}
        for k in dropped[:3]:                              # the last dropped parameter stays missing
            src = rng.choice(pool[len(k)])
            analogs[k if rng.random() < 0.5 else k[::-1]] = (src if rng.random() < 0.5 else src[::-1],
                                                            rng.choice([0.6, 12.0, 49.0, 51.0, 120.0]))
        r = assign_charges(types, bonds, [0] * len(types), IncrementTable(sub_raw), analogs)
        pq, ppen, pmissing = proto.assign(types, bonds, [0] * len(types), sub_proto, analogs)
        assert r.charges == pq and r.penalties == ppen and r.missing == pmissing, res.name


# ---- unit tests on a synthetic table -------------------------------------------------------------------------------

def _table():
    return IncrementTable({2: {('A', 'B'): (0.1,), ('B', 'C'): None},
                           3: {('A', 'B', 'C'): (0.01, 0.02)},
                           4: {('A', 'B', 'C', 'D'): (0.001, 0.002, 0.003)}})


def test_bond_increment_moves_charge_from_first_to_second_atom():
    r = assign_charges(['A', 'B'], [(0, 1)], [0, 0], _table())
    assert r.charges == pytest.approx([-0.1, 0.1]) and r.penalties == [0.0, 0.0] and not r.missing


def test_reversed_instance_is_reversed():
    r = assign_charges(['B', 'A'], [(0, 1)], [0, 1], _table())
    assert r.charges == pytest.approx([0.1, 0.9])


def test_na_increments_contribute_nothing_and_are_not_missing():
    r = assign_charges(['B', 'C'], [(0, 1)], [0, 0], _table())
    assert r.charges == [0.0, 0.0] and not r.missing


def test_angle_and_dihedral_increments_follow_bonds():
    t = IncrementTable({2: {('A', 'B'): (0.0,), ('B', 'C'): (0.0,), ('C', 'D'): (0.0,)},
                        3: {('A', 'B', 'C'): (0.01, 0.02), ('B', 'C', 'D'): (0.0, 0.0)},
                        4: {('A', 'B', 'C', 'D'): (0.001, 0.002, 0.003)}})
    r = assign_charges(['A', 'B', 'C', 'D'], [(0, 1), (1, 2), (2, 3)], [0] * 4, t)
    assert r.charges == pytest.approx([-0.011, 0.01 - 0.02 + 0.001 - 0.002, 0.02 + 0.002 - 0.003, 0.003])


def test_missing_parameters_are_reported():
    r = assign_charges(['A', 'X'], [(0, 1)], [0, 0], _table())
    assert r.missing == [('A', 'X')] and r.charges == [0.0, 0.0]


def test_analog_supplies_increment_and_penalty_in_either_orientation():
    for analogs in ({('X', 'A'): Analog(('B', 'A'), 2.0)}, {('A', 'X'): (('A', 'B'), 2.0)}):
        r = assign_charges(['X', 'A'], [(0, 1)], [0, 0], _table(), analogs)
        assert r.charges == pytest.approx([0.1, -0.1]) and not r.missing
        assert r.penalties == pytest.approx([atom_penalty([(0.1, 2.0)]), atom_penalty([(-0.1, 2.0)])])
        assert list(r.analogs_used.values())[0].penalty == 2.0


def test_high_penalty_dihedral_analog_is_ignored():
    t = _table()
    for pen, expect_q, expect_p in ((50.0, 0.001, 50.0), (50.5, 0.0, 50.0)):
        r = assign_charges(['A', 'B', 'C', 'E'], [(0, 1), (1, 2), (2, 3)], [0] * 4,
                           IncrementTable({**t._table, 3: {**t._table[3], ('B', 'C', 'E'): (0.0, 0.0)},
                                           2: {**t._table[2], ('C', 'E'): (0.0,)}}),
                           {('A', 'B', 'C', 'E'): (('A', 'B', 'C', 'D'), pen)})
        assert r.charges[3] == pytest.approx(0.003 if expect_q else 0.0)
        assert r.penalties[3] == pytest.approx(atom_penalty([(r.charges[3], expect_p)]))


def test_atom_penalty_formula():
    assert atom_penalty([]) == 0.0
    assert atom_penalty([(0.0, 3.0)]) == pytest.approx(math.sqrt((0.05 ** 6) ** (1 / 3) * 9))
    assert atom_penalty([(0.1, 1.0), (-0.2, 2.0)]) == pytest.approx(
        math.sqrt((0.1 + 0.05 ** 6) ** (1 / 3) + 4 * (0.2 + 0.05 ** 6) ** (1 / 3)))


def test_typed_molecule_validation():
    tm = TypedMolecule(['C1', 'O1'], ['CG2O3', 'OG2D2'], [(0, 1)], [-1, 0])
    assert tm.net_charge == -1
    assert TypedMolecule(['A'], ['HGA1'], []).formal_charges == [0]
    with pytest.raises(ValueError):
        TypedMolecule(['A', 'B'], ['X'], [])
    with pytest.raises(ValueError):
        TypedMolecule(['A', 'B'], ['X', 'Y'], [(0, 2)])


def test_prm_reader_rejects_wrong_increment_count(tmp_path):
    p = tmp_path / 'bad.prm'
    p.write_text('BONDS\nA B 100.0 1.5 !<% 0.1 0.2 %>!\nEND\n')
    with pytest.raises(ValueError):
        read_increments(p)

"""Typing rule language (SPEC 5.5), the P-I Fig. 5 case studies, and the production rule file (SPEC 5.6-5.7)."""
import os

import pytest
from rdkit import Chem

from open_cgenff.molecule import InputMolecule
from open_cgenff.perception import perceive
from open_cgenff.typing import (RULES_FILE, RuleSet, RuleSyntaxError, TypingError, default_rules, parse,
                                type_molecule, type_perception)

FIXTURES = os.path.join(os.path.dirname(__file__), 'fixtures')


def perceived(smiles):
    return perceive(InputMolecule.from_rdkit(Chem.AddHs(Chem.MolFromSmiles(smiles))))


def run(rules_text, smiles):
    return type_perception(parse(rules_text), perceived(smiles))


def types(rules_text, smiles, heavy_only=True):
    p = perceived(smiles)
    res = type_perception(parse(rules_text), p)
    return [t for t, e in zip(res.types, p.molecule.elements) if not heavy_only or e != 'H']


def one_rule(conditions, yes='Y', no='N'):
    """Rules that type an atom `yes` when conditions hold (H atoms always 'H')."""
    return f"cat main\ntyp H : el H\ntyp {yes} : {conditions}\ntyp {no} :\n"


# ---- parser -----------------------------------------------------------------------------------------------------

def test_parser_categories_macros_and_comments():
    rs = parse('def two ( ne ( el H ) ( el H ) )\ncat main\nsub X : @two   # comment\ntyp A :\ncat X\n'
               'typ B : warn "a # in a message"\nend\n')
    assert list(rs.categories) == ['main', 'X']
    assert rs.categories['main'][0].conditions == [('ne', [[('el', 'H')], [('el', 'H')]])]
    assert rs.categories['X'][0].warn == ['a # in a message']
    assert rs.types() == ['A', 'B']


def test_parser_multiline_rule():
    rs = parse('cat main\ntyp A : ne ( el C\n   nh 3 )\ntyp B :\n')
    assert rs.categories['main'][0].conditions == [('ne', [[('el', 'C'), ('nh', 3)]])]


@pytest.mark.parametrize('text', [
    'cat main\ntyp A : foo\n',                       # unknown keyword
    'cat main\ntyp A : ne ( el C\n',                 # unclosed parenthesis
    'cat main\nsub X :\n',                           # sub to an undefined category
    'cat other\ntyp A :\n',                          # no cat main
    'typ A :\n',                                     # rule outside a category
    'cat main\ntyp A : ne ( charge 1 )\n',           # action inside a group
    'cat main\ntyp A : @nope\n',                     # undefined macro
])
def test_parser_errors(text):
    with pytest.raises(RuleSyntaxError):
        parse(text)


# ---- one test per keyword ---------------------------------------------------------------------------------------

def test_el_elha_elos():
    assert types(one_rule('el O'), 'CO') == ['N', 'Y']
    assert types(one_rule('elha'), 'C(F)(Cl)(Br)I') == ['N', 'Y', 'Y', 'Y', 'Y']
    assert types(one_rule('elos'), 'OCCS') == ['Y', 'N', 'N', 'Y']


def test_nb_counts_bond_orders_including_h():
    assert types(one_rule('nb 4'), 'C=CC#N') == ['Y', 'Y', 'Y', 'N']
    assert types(one_rule('nb 4'), 'C[NH3+]') == ['Y', 'Y']           # ammonium N: 4 single bonds


def test_rings_and_ring_classes():
    assert types(one_rule('rings 2'), 'c1ccc2ccccc2c1')[3] == 'Y'      # fusion atom
    assert types(one_rule('arom 6'), 'c1ccccc1C') == ['Y'] * 6 + ['N']
    assert types(one_rule('ring3 6'), 'C1CCCCC1')[0] == 'Y'           # all-sp3 ring
    assert types(one_rule('ring2 5'), 'O=C1C=CC(=O)N1')[1] == 'Y'      # maleimide: all-sp2, not aromatic
    assert types(one_rule('ring23 6'), 'C1=CCCCC1')[0] == 'Y'          # mixed ring
    assert types(one_rule('ring 5'), 'C1CCCC1C') == ['Y'] * 5 + ['N']
    assert types(one_rule('arom'), 'c1ccoc1')[0] == 'Y'               # any size


def test_matched_ring_is_consumed():
    rule = one_rule('arom 6 arom 6')
    t = types(rule, 'c1ccc2ccccc2c1')
    assert t[3] == 'Y' and t[8] == 'Y' and t[0] == 'N'


def test_ne_greedy_order_dependence():
    # N(C0 methyl, C2 methylene): the first group takes the lowest-index C; no backtracking (P-I)
    assert types(one_rule('el N ne ( el C ) ( el C nh 3 )'), 'CN(C)C')[1] == 'Y'
    assert types(one_rule('el N ne ( el C ) ( el C nh 3 )'), 'CNCC')[1] == 'N'
    assert types(one_rule('el N ne ( el C nh 3 ) ( el C )'), 'CNCC')[1] == 'Y'


def test_ne_does_not_go_back_and_self_closes_rings():
    assert types(one_rule('el C ne ( ne ( el C ) )'), 'CC') == ['N', 'N']          # parent not revisited
    assert types(one_rule('el C ne ( ne ( ne ( self ) ) )'), 'C1CC1') == ['Y'] * 3  # 3-ring closure
    assert types(one_rule('el C ne ( ne ( ne ( self ) ) )'), 'C1CCC1') == ['N'] * 4


def test_bo_and_inring():
    assert types(one_rule('ne ( bo 2 el O )'), 'CC=O') == ['N', 'Y', 'N']
    assert types(one_rule('bo 1'), 'CC') == ['N', 'N']                             # no bond at top level
    assert types(one_rule('ne ( el C inring )'), 'C1CC1C') == ['Y', 'Y', 'Y', 'N']
    assert types(one_rule('ne ( el C ! ( inring ) )'), 'C1CC1C') == ['N', 'N', 'Y', 'Y']
    assert types(one_rule('ne ( inring arom 6 )'), 'c1ccccc1C1CC1')[6] == 'N'     # extension: class/size
    assert types(one_rule('ne ( inring 3 )'), 'c1ccccc1C1CC1')[6] == 'Y'


def test_not_and_or():
    assert types(one_rule('! ( el C )'), 'CO') == ['N', 'Y']
    assert types(one_rule('or ( el O ) ( el N )'), 'CON') == ['N', 'Y', 'Y']


def test_extensions_nh_deg_fc_type_rmember_fusion():
    assert types(one_rule('nh 2'), 'CCC') == ['N', 'Y', 'N']
    assert types(one_rule('deg 1'), 'C#N') == ['N', 'Y']
    assert types(one_rule('fc -1'), 'CC(=O)[O-]')[3] == 'Y'
    assert types('cat main\ntyp H : el H ne ( type C2 )\ntyp X : el H\ntyp C2 : nh 2\ntyp C :\n', 'CC=C',
                 heavy_only=False)[-2:] == ['H', 'H']
    assert types(one_rule('rmember arom 6 ( el N )'), 'c1ccncc1C')[:6] == ['Y'] * 6
    assert types(one_rule('fusion'), 'C1CCC2CCCCC2C1')[3] == 'Y'
    assert types(one_rule('fusion'), 'C1CCC2(CC1)CCCC2')[3] == 'N'                 # spiro atom


def test_actions_charge_impr_warn_err_and_trace():
    res = run('cat main\nsub X : el N charge 1 warn "w1"\ntyp H : el H\ntyp C : impr\ncat X\ntyp N4 : err "e1"\n',
              'C=[NH2+]')
    assert res.types[:2] == ['C', 'N4'] and res.charges[:2] == [0, 1]
    assert res.impropers == [(0, 1, 2, 3)]
    assert any('w1' in w for w in res.warnings) and any('e1' in e for e in res.errors)
    assert [r.target for r in res.atoms[1].trace] == ['X', 'N4']


def test_no_matching_rule_and_type_cycles():
    res = run('cat main\ntyp A : el O\n', 'C')
    assert res.types[0] == 'UNK' and res.errors
    with pytest.raises(TypingError):
        run('cat main\ntyp A : ne ( type B )\ntyp B :\n', 'CC')


def test_altnum_numbering():
    rules = 'cat main\ntyp H : el H\ntyp CG2DC? : ne ( bo 2 ) altnum\ntyp C :\n'
    assert types(rules, 'C=CC=CC=C') == ['CG2DC1', 'CG2DC1', 'CG2DC2', 'CG2DC2', 'CG2DC1', 'CG2DC1']


# ---- P-I Fig. 5 case studies (SPEC 5.5 acceptance) -------------------------------------------------------------

@pytest.fixture(scope='module')
def fig5():
    return RuleSet.from_file(os.path.join(FIXTURES, 'pi_fig5.rules'))


def test_fig5_pyridinium(fig5):
    p = perceived('c1cc[nH+]cc1')
    res = type_perception(fig5, p)
    n = p.molecule.elements.index('N')
    assert res.types[n] == 'NG2R61' and res.charges[n] == 1


def test_fig5_pyridine_and_pyrimidine(fig5):
    # NG2R62 needs another ring N ortho/meta; the rule only works because nested ne does not revisit the parent
    for smi, want in (('c1ccncc1', 'NG2R60'), ('c1cncnc1', 'NG2R62'), ('c1ccnnc1', 'NG2R62'), ('c1cnccn1', 'NG2R60')):
        p = perceived(smi)
        res = type_perception(fig5, p)
        assert {t for t, e in zip(res.types, p.molecule.elements) if e == 'N'} == {want}, smi


def test_fig5_pentadiene_altnum(fig5):
    t = types_with(fig5, 'C=CC=CC')
    assert t[:4] == ['CG2DC3', 'CG2DC1', 'CG2DC2', 'CG2DC2']


def types_with(rules, smiles):
    p = perceived(smiles)
    return [t for t, e in zip(type_perception(rules, p).types, p.molecule.elements) if e != 'H']


def test_fig5_iminium_charge_and_impropers(fig5):
    p = perceived('CC=[NH+]C')
    res = type_perception(fig5, p)
    assert res.charges[1] == 1 and res.impropers and res.impropers[0][0] == 1


# ---- production rule file -----------------------------------------------------------------------------------------

def test_production_rules_parse_and_use_known_types(toppar):
    rules = default_rules()
    known = set()
    for line in open(os.path.join(toppar, 'top_all36_cgenff.rtf')):
        if line.startswith('MASS'):
            known.add(line.split()[2])
    for t in rules.types():
        variants = [t] if '?' not in t else [t.replace('?', '1'), t.replace('?', '2')]
        assert all(v in known for v in variants), t
    assert RULES_FILE.name == 'cgenff_v5.rules'


@pytest.mark.parametrize('smiles,expected', [
    ('CC(=O)[O-]', ['CG331', 'CG2O3', 'OG2D2', 'OG2D2']),
    ('C[NH3+]', ['CG334', 'NG3P3']),
    ('c1ccccc1Cl', ['CG2R61'] * 5 + ['CG2R61', 'CLGR1']),
    ('NC(=[NH2+])N', ['NG2P1', 'CG2N1', 'NG2P1', 'NG2P1']),
    ('C=CC=C', ['CG2DC3', 'CG2DC1', 'CG2DC2', 'CG2DC3']),
    ('c1ccc2c(c1)cccc2', ['CG2R61'] * 10),                         # 6/6 fusion is not CG2RC0
])
def test_production_typer_examples(smiles, expected):
    tm = type_molecule(InputMolecule.from_rdkit(Chem.AddHs(Chem.MolFromSmiles(smiles))))
    heavy = [t for t in tm.types if not t.startswith('H')]
    assert heavy == expected


def test_production_formal_charges():
    for smi, q in (('CC(=O)[O-]', -1), ('C[NH3+]', 1), ('NC(=[NH2+])N', 1), ('c1cc[nH+]cc1', 1),
                   ('COP(=O)([O-])[O-]', -2), ('CS(=O)(=O)[O-]', -1), ('C[O-]', -1)):
        tm = type_molecule(InputMolecule.from_rdkit(Chem.AddHs(Chem.MolFromSmiles(smi))))
        assert tm.net_charge == q, (smi, tm.types, tm.formal_charges)

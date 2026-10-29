"""Benchmark B1: typing reference from the CGenFF 5.0 topology file (SPEC 7.1, FINDINGS 3.1).

For each RESI with a SMILES comment, the RTF atoms are mapped 1:1 onto an RDKit molecule (explicit H, bond orders)
by graph isomorphism; only exact bond-graph matches are kept. SMILES whose bond orders were lost (radical atoms) are
repaired from connectivity + net charge. Split membership is read from the fixed manifest
(benchmarks/b1_rtf_zinc.split.json); `shuffle_split` reproduces how that manifest was made (seed 42, 80/20).
"""
import json
import random

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds

from ..io.rtf import read_rtf
from ..paths import REPO_ROOT
from .common import TypingCase, element_of_type

RDLogger.DisableLog("rdApp.*")
MANIFEST = REPO_ROOT / "benchmarks" / "b1_rtf_zinc.split.json"


def _repair_bond_orders(mol, net_charge):
    m = Chem.RWMol(mol)
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetNumRadicalElectrons(0); a.SetNoImplicit(True); a.SetIsAromatic(False)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    try:
        rdDetermineBonds.DetermineBondOrders(m, charge=int(round(net_charge)), embedChiral=False)
        m = m.GetMol()
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return None if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()) else m


def _map_residue(res):
    """RDKit molecule + reference types (RDKit order) for one residue, or (None, reason)."""
    mol = Chem.MolFromSmiles(res.smiles)
    if mol is None:
        return None, "bad SMILES"
    mol = Chem.AddHs(mol)
    types = res.atom_types()
    names = [a.name for a in res.atoms if element_of_type(a.type)]
    if len(names) != mol.GetNumAtoms():
        return None, "atom count"
    idx = {n: i for i, n in enumerate(names)}
    query = Chem.RWMol()
    for n in names:
        query.AddAtom(Chem.AtomFromSmarts(f"[#{Chem.Atom(element_of_type(types[n])).GetAtomicNum()}]"))
    seen = set()
    for a, b in res.bonds:
        if a in idx and b in idx and frozenset((a, b)) not in seen:
            seen.add(frozenset((a, b)))
            query.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
            query.ReplaceBond(query.GetBondBetweenAtoms(idx[a], idx[b]).GetIdx(), Chem.BondFromSmarts("~"))
    if len(seen) != mol.GetNumBonds():
        return None, "bond count"
    match = mol.GetSubstructMatch(query.GetMol())
    if len(match) != mol.GetNumAtoms():
        return None, "bond graph"
    repaired = False
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        mol = _repair_bond_orders(mol, res.net_charge)
        if mol is None:
            return None, "unrepairable radical"
        repaired = True
    ref, charges = [None] * mol.GetNumAtoms(), [None] * mol.GetNumAtoms()
    q = {a.name: a.charge for a in res.atoms}
    for qi, mi in enumerate(match):
        ref[mi] = types[names[qi]]
        charges[mi] = q[names[qi]]
    impropers = [tuple(match[idx[n]] for n in imp) for imp in res.impropers if all(n in idx for n in imp)]
    return (mol, ref, charges, impropers, repaired), None


def build(rtf_path, manifest=MANIFEST):
    """Return (cases, stats). With manifest=None, every case gets split 'unassigned' (see shuffle_split)."""
    split_of = {}
    if manifest is not None:
        for split, names in json.loads(open(manifest).read()).items():
            split_of.update({n: split for n in names})
    cases, skipped = [], {}
    n_repaired = 0
    for res in read_rtf(rtf_path):
        if res.kind != "RESI" or not res.smiles:
            continue
        out, reason = _map_residue(res)
        if out is None:
            skipped[reason] = skipped.get(reason, 0) + 1
            continue
        mol, ref, charges, impropers, repaired = out
        n_repaired += repaired
        cases.append(TypingCase(res.name, mol, ref, split_of.get(res.name, "unassigned"), res.smiles,
                                {"origin": "zinc" if res.name.startswith("230") else "legacy",
                                 "lonepair": res.has_lonepair, "charges": charges,
                                 "impropers": impropers}))
    if manifest is not None:
        built = {c.name for c in cases}
        if built != set(split_of) or any(c.split == "unassigned" for c in cases):
            raise ValueError(f"B1 build does not match manifest {manifest}: "
                             f"{len(built - set(split_of))} extra, {len(set(split_of) - built)} missing")
    stats = {"molecules": len(cases), "atoms": sum(len(c.types) for c in cases), "repaired_radicals": n_repaired,
             "skipped": skipped}
    return cases, stats


def shuffle_split(cases, seed=42):
    """The procedure that produced the manifest: shuffle in file order, first 20 % = test. Returns {split: [names]}."""
    order = list(cases)
    random.Random(seed).shuffle(order)
    n_test = len(order) // 5
    return {"train": sorted(c.name for c in order[n_test:]), "test": sorted(c.name for c in order[:n_test])}

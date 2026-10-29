"""Benchmark B2: typing reference from the original CGenFF model compounds (RTF residues without SMILES; SPEC 7.1).

Bond orders come from connectivity + the residue's net charge (RDKit DetermineBondOrders); elements from the atom
types; lone-pair sites (LPH) are dropped. The IC tables are not needed for bond orders. Residues whose bond orders
cannot be determined are skipped. Split membership is read from benchmarks/b2_rtf_legacy.split.json (made once with
the same procedure as B1: seed 42, 80/20; see b1.shuffle_split).
"""
import json

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds

from ..io.rtf import read_rtf
from ..paths import REPO_ROOT
from .common import TypingCase, element_of_type

RDLogger.DisableLog("rdApp.*")
MANIFEST = REPO_ROOT / "benchmarks" / "b2_rtf_legacy.split.json"


def _neutralize(m):
    """Remove charge separation that DetermineBondOrders sometimes produces (e.g. a thioester as [S+]=C-[O-]):
    while some +1 atom starts an alternating double/single path that ends with a single bond into a -1 atom, flip
    the bond orders along it and clear both charges. Works on a kekulized RWMol; returns the number of flips."""
    flips = 0
    while True:
        path = None
        for a in m.GetAtoms():
            if a.GetFormalCharge() != 1:
                continue
            # depth-first over alternating paths: (atom, bond list, next order wanted)
            stack = [(a.GetIdx(), [], Chem.BondType.DOUBLE)]
            while stack and path is None:
                i, bonds, want = stack.pop()
                for b in m.GetAtomWithIdx(i).GetBonds():
                    j = b.GetOtherAtomIdx(i)
                    if b.GetBondType() != want or b.GetIdx() in bonds or j == a.GetIdx() or len(bonds) > 12:
                        continue
                    if want == Chem.BondType.SINGLE and m.GetAtomWithIdx(j).GetFormalCharge() == -1:
                        path = (a.GetIdx(), j, bonds + [b.GetIdx()])
                        break
                    nxt = Chem.BondType.SINGLE if want == Chem.BondType.DOUBLE else Chem.BondType.DOUBLE
                    stack.append((j, bonds + [b.GetIdx()], nxt))
            if path:
                break
        if path is None:
            return flips
        start, end, bonds = path
        for k in bonds:
            b = m.GetBondWithIdx(k)
            b.SetBondType(Chem.BondType.SINGLE if b.GetBondType() == Chem.BondType.DOUBLE else Chem.BondType.DOUBLE)
        for k in (start, end):
            m.GetAtomWithIdx(k).SetFormalCharge(0)
        flips += 1


def _build_mol(res):
    """(RDKit molecule, reference types), or (None, reason)."""
    atoms = [a for a in res.atoms if element_of_type(a.type)]
    if len(atoms) < 2:
        return None, "fewer than 2 atoms"
    idx = {a.name: i for i, a in enumerate(atoms)}
    m = Chem.RWMol()
    for a in atoms:
        at = Chem.Atom(element_of_type(a.type))
        at.SetNoImplicit(True)
        m.AddAtom(at)
    for a, b in res.bonds:
        if a in idx and b in idx and a != b and m.GetBondBetweenAtoms(idx[a], idx[b]) is None:
            m.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
    if len(Chem.GetMolFrags(m)) != 1:
        return None, "disconnected"
    try:
        rdDetermineBonds.DetermineBondOrders(m, charge=int(round(res.net_charge)), embedChiral=False)
        Chem.Kekulize(m, clearAromaticFlags=True)
        _neutralize(m)
        mol = m.GetMol()
        Chem.SanitizeMol(mol)
    except Exception:
        return None, "bond orders"
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None, "radical"
    return (mol, [a.type for a in atoms], [a.charge for a in atoms]), None


def build(rtf_path, manifest=MANIFEST):
    """Return (cases, stats). With manifest=None every case gets split 'unassigned' (see b1.shuffle_split)."""
    split_of = {}
    if manifest is not None:
        for split, names in json.loads(open(manifest).read()).items():
            split_of.update({n: split for n in names})
    cases, skipped = [], {}
    for res in read_rtf(rtf_path):
        if res.kind != "RESI" or res.smiles:
            continue
        out, reason = _build_mol(res)
        if out is None:
            skipped[reason] = skipped.get(reason, 0) + 1
            continue
        mol, ref, charges = out
        cases.append(TypingCase(res.name, mol, ref, split_of.get(res.name, "unassigned"),
                                Chem.MolToSmiles(Chem.RemoveHs(mol)),
                                {"origin": "legacy", "lonepair": res.has_lonepair, "net_charge": res.net_charge,
                                 "charges": charges}))
    if manifest is not None:
        built = {c.name for c in cases}
        if built != set(split_of):
            raise ValueError(f"B2 build does not match manifest {manifest}: "
                             f"{len(built - set(split_of))} extra, {len(set(split_of) - built)} missing")
    stats = {"molecules": len(cases), "atoms": sum(len(c.types) for c in cases), "skipped": skipped}
    return cases, stats

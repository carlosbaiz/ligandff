"""Benchmark B3: CHARMM-GUI GROMACS topologies with reference CGenFF types (CGenFF <= 4.6 era; SPEC 7.1).

RDKit molecules are built from the topologies (elements from masses, bond orders from connectivity + net charge).
Only topologies are used: <A>.itp, or cgen_<A>.itp where the two files are swapped (mbn, au); parameter files are not
read and not committed. Known data issues: benchmarks/charmm_gui_cpa/NOTES.md.
"""
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds

from ..paths import REPO_ROOT
from .common import TypingCase

RDLogger.DisableLog("rdApp.*")
RAW_DIR = REPO_ROOT / "benchmarks" / "charmm_gui_cpa" / "raw"
MASSES = [(1.008, "H"), (10.811, "B"), (12.011, "C"), (14.007, "N"), (15.9994, "O"), (18.998, "F"), (30.974, "P"),
          (32.06, "S"), (35.45, "Cl"), (79.904, "Br"), (126.904, "I"), (196.97, "Au"), (22.99, "Na"), (78.96, "Se")]


def _element(mass):
    return min(MASSES, key=lambda x: abs(x[0] - mass))[1] if mass > 0.5 else None


def _parse_itp(path):
    section, name, atoms, bonds = None, None, [], []
    for line in open(path, errors="ignore"):
        s = line.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            section = s.strip("[] ").lower()
            continue
        w = s.split()
        if section == "moleculetype" and name is None:
            name = w[0]
        elif section == "atoms" and len(w) >= 8:
            atoms.append({"type": w[1], "name": w[4], "q": float(w[6]), "mass": float(w[7])})
        elif section == "bonds" and len(w) >= 2:
            bonds.append((int(w[0]) - 1, int(w[1]) - 1))
    return name, atoms, bonds


def _build_mol(atoms, bonds):
    keep = [i for i, a in enumerate(atoms) if _element(a["mass"])]      # drop virtual sites / lone pairs (mass 0)
    idx = {i: k for k, i in enumerate(keep)}
    m = Chem.RWMol()
    for i in keep:
        a = Chem.Atom(_element(atoms[i]["mass"]))
        a.SetNoImplicit(True)
        m.AddAtom(a)
    for i, j in bonds:
        if i in idx and j in idx and m.GetBondBetweenAtoms(idx[i], idx[j]) is None:
            m.AddBond(idx[i], idx[j], Chem.BondType.SINGLE)
    net = round(sum(a["q"] for a in atoms))
    if net == 0 and any(atoms[i]["type"] == "SG302" for i in keep):
        net = -1                                                        # MBN: thiolate with custom neutral charges
    rdDetermineBonds.DetermineBondOrders(m, charge=int(net), embedChiral=False)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return mol, keep, net


def _topology(raw_dir, species):
    for fn in (raw_dir / f"{species}.itp", raw_dir / f"cgen_{species}.itp"):
        if fn.exists():
            name, atoms, bonds = _parse_itp(fn)
            if atoms:
                return fn.name, name, atoms, bonds
    return None, None, [], []


def build(raw_dir=RAW_DIR):
    """Return (cases, excluded): excluded = [(species, reason)]."""
    raw_dir = Path(raw_dir)
    species = sorted({p.stem.removeprefix("cgen_") for p in raw_dir.glob("*.itp")})
    cases, excluded = [], []
    for sp in species:
        fname, name, atoms, bonds = _topology(raw_dir, sp)
        if not atoms:
            excluded.append((sp, "no atoms section"))
            continue
        try:
            mol, keep, net = _build_mol(atoms, bonds)
        except Exception as e:
            excluded.append((sp, f"could not build molecule ({type(e).__name__})"))
            continue
        if mol.GetNumAtoms() < 2 or (mol.GetNumAtoms() <= 3 and all(a.GetSymbol() in ("H", "O") for a in mol.GetAtoms())):
            excluded.append((sp, "water / single atom (not a CGenFF ligand)"))
            continue
        ref = [atoms[i]["type"] for i in keep]
        if any(t[:2] in ("CC", "HC", "OC") and not t.startswith(("CG", "HG", "OG", "CL")) for t in ref):
            excluded.append((sp, "CHARMM36 carbohydrate types, not CGenFF"))
            continue
        cases.append(TypingCase(name, mol, ref, "dev", Chem.MolToSmiles(Chem.RemoveHs(mol)),
                                {"file": fname, "charges": [atoms[i]["q"] for i in keep], "net_charge": net}))
    return cases, excluded

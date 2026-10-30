"""Tripos mol2 reader (SPEC 5.1). Keeps input atom order; 'ar' bonds are left for resonance resolution.

Element: the SYBYL type before the dot ('C.ar' -> C, 'N.4' -> N; 'H1' -> H), else the atom name's letters.
Bond types: '1', '2', '3', 'ar'; 'am' (amide) is a single bond; 'du'/'un'/'nc' are rejected. Formal charges are
not stored in mol2 (the charge column holds partial charges), so they are left to perception.
"""
import re

from ..molecule import AROMATIC, InputMolecule

BOND_ORDERS = {"1": 1, "2": 2, "3": 3, "ar": AROMATIC, "am": 1}
_ELEMENTS = {"H", "B", "C", "N", "O", "F", "P", "S", "CL", "BR", "I", "SE", "AL", "SI", "NA", "K", "LI", "MG", "CA",
             "ZN", "FE", "CU"}


def _element(sybyl, name):
    token = re.sub(r"\d+$", "", sybyl.split(".")[0]).upper()
    if token in _ELEMENTS:
        return token.capitalize()
    letters = re.sub(r"[^A-Za-z]", "", name).upper()
    for k in (letters[:2], letters[:1]):
        if k in _ELEMENTS and (len(k) == 1 or k in ("CL", "BR")):
            return k.capitalize()
    raise ValueError(f"cannot determine element of atom {name!r} (type {sybyl!r})")


def read_mol2(path):
    """Return the molecules (InputMolecule) of a mol2 file."""
    out, cur, section = [], None, None
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>"):
                section = line.strip()[9:]
                if section == "MOLECULE":
                    cur = {"name": None, "atoms": [], "bonds": []}
                    out.append(cur)
                continue
            w = line.split()
            if not w or line.startswith("#") or cur is None:
                continue
            if section == "MOLECULE" and cur["name"] is None:
                cur["name"] = line.strip()
            elif section == "ATOM" and len(w) >= 6:
                cur["atoms"].append((w[0], w[1], float(w[2]), float(w[3]), float(w[4]), w[5]))
            elif section == "BOND" and len(w) >= 4:
                if w[3] not in BOND_ORDERS:
                    raise ValueError(f"{path}: unsupported bond type {w[3]!r}")
                cur["bonds"].append((w[1], w[2], BOND_ORDERS[w[3]]))
    mols = []
    for m in out:
        index = {a[0]: k for k, a in enumerate(m["atoms"])}
        mols.append(InputMolecule(m["name"] or "", [a[1] for a in m["atoms"]], [_element(a[5], a[1]) for a in m["atoms"]],
                                  [(index[a], index[b], o) for a, b, o in m["bonds"]], None,
                                  [a[2:5] for a in m["atoms"]]))
    return mols

"""Halogen lone pairs (P-X; SPEC 5.12).

Every aromatic Cl/Br/I (types CLGR1, BRGR1, IGR1) gets a massless, LJ-free virtual site LPH, colinear with the C-X
bond, at 1.64 / 1.89 / 2.20 A from X, charge +0.05 e. The halogen keeps the molecule neutral: its increment charge
minus 0.05 (FINDINGS 2.3: best simple rule on the training residues, max error 0.04 e). The RTF has two aliphatic
halogens with lone pairs (SM200, SM211); P-X does not, and neither does this module.
"""
from dataclasses import dataclass

LP_TYPE = "LPH"
LP_CHARGE = 0.05
LP_DISTANCE = {"CLGR1": 1.64, "BRGR1": 1.89, "IGR1": 2.20}


@dataclass(frozen=True)
class LonePair:
    name: str
    host: int                 # halogen atom index
    base: int                 # carbon bonded to the halogen
    distance: float
    charge: float = LP_CHARGE
    type: str = LP_TYPE


def lone_pairs(names, types, bonds):
    """Lone pairs for a typed molecule, named LP1, LP2, ... in atom order (skipping names already in use)."""
    neighbours = {}
    for i, j in bonds:
        neighbours.setdefault(i, []).append(j)
        neighbours.setdefault(j, []).append(i)
    used = set(names)
    out, k = [], 0
    for i, t in enumerate(types):
        if t not in LP_DISTANCE or len(neighbours.get(i, ())) != 1:
            continue
        k += 1
        while f"LP{k}" in used:
            k += 1
        name = f"LP{k}"
        used.add(name)
        out.append(LonePair(name, i, neighbours[i][0], LP_DISTANCE[t]))
    return out


def apply_lone_pair_charges(charges, lps):
    """Halogen charges after moving +0.05 e onto each lone pair (returns a new list)."""
    q = list(charges)
    for lp in lps:
        q[lp.host] -= lp.charge
    return q

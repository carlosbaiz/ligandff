"""Reader for CHARMM topology (RTF) files: residues with atoms, bonds and the SMILES comments of CGenFF 5.0 entries."""
from dataclasses import dataclass, field


@dataclass
class RtfAtom:
    name: str
    type: str
    charge: float


@dataclass
class Residue:
    name: str
    kind: str                                   # 'RESI' or 'PRES'
    net_charge: float
    atoms: list = field(default_factory=list)   # [RtfAtom] in file order
    bonds: list = field(default_factory=list)   # [(name, name)] from BOND/DOUBLE/TRIPLE, duplicates kept
    smiles: str = None                          # first 'SMILES:' comment after the RESI line
    has_lonepair: bool = False
    has_external_bond: bool = False             # a bond to a +/- atom of a neighbouring residue

    def atom_types(self):
        return {a.name: a.type for a in self.atoms}


def read_rtf(path):
    """Return the residues (RESI and PRES) of an RTF file, up to the first END card."""
    residues, cur = [], None
    with open(path, errors="ignore") as fh:
        for line in fh:
            words = line.split("!")[0].split()
            if words and words[0] in ("RESI", "PRES"):
                cur = Residue(words[1], words[0], float(words[2]))
                residues.append(cur)
                continue
            if cur is None:
                continue
            if "SMILES:" in line and cur.smiles is None:
                cur.smiles = line.split("SMILES:")[1].split()[0]
            if not words:
                continue
            if words[0] == "ATOM":
                cur.atoms.append(RtfAtom(words[1], words[2], float(words[3])))
            elif words[0] in ("BOND", "DOUBLE", "TRIPLE"):
                pairs = list(zip(words[1::2], words[2::2]))
                cur.bonds += pairs
                cur.has_external_bond |= any(a[0] in "+-" or b[0] in "+-" for a, b in pairs)
            elif words[0] in ("LONEPAIR", "LONE"):
                cur.has_lonepair = True
            elif words[0] == "END":
                break
    return residues

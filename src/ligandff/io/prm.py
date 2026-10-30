"""Reader for the charge-increment field of CGenFF 5.0 parameter files (par_all36_cgenff_increments.prm).

Each bond/angle/dihedral line may carry `!<% d1 [d2 [d3]] %>!`: one increment per bond of the parameter, in the
parameter's atom order; `N/A` = no increment (the term does not contribute to charges). Multi-term dihedrals repeat
the same increments on every line; the first line counts.
"""
import re
from dataclasses import dataclass, field

SECTIONS = {"BONDS": 2, "ANGLES": 3, "DIHEDRALS": 4}
OTHER_SECTIONS = ("IMPROPER", "IMPROPERS", "NONBONDED", "NBFIX", "HBOND", "END", "CMAP")
_FIELD = re.compile(r"<%(.*?)%>")


def read_increments(path):
    """Return {n_atoms: {types: tuple of increments, or None for N/A}} for n_atoms = 2, 3, 4."""
    table = {n: {} for n in SECTIONS.values()}
    section = None
    with open(path, errors="ignore") as fh:
        for line in fh:
            words = line.split("!")[0].split()
            if words and words[0].upper() in (*SECTIONS, *OTHER_SECTIONS):
                section = words[0].upper()
                continue
            if section not in SECTIONS or not words:
                continue
            m = _FIELD.search(line)
            if not m:
                continue
            n = SECTIONS[section]
            key = tuple(words[:n])
            if key in table[n]:
                continue
            body = m.group(1)
            if "N/A" in body:
                table[n][key] = None
                continue
            values = tuple(float(x) for x in body.split())
            if len(values) != n - 1:
                raise ValueError(f"{path}: expected {n - 1} increments for {' '.join(key)}, got {len(values)}")
            table[n][key] = values
    return table


# ---- full parameter file (SPEC 5.8) --------------------------------------------------------------------------------

BONDED = {"BONDS": "bond", "ANGLES": "angle", "DIHEDRALS": "dihedral", "IMPROPER": "improper", "IMPROPERS": "improper"}
N_ATOMS = {"bond": 2, "angle": 3, "dihedral": 4, "improper": 4}
_ANNOTATION = re.compile(r"from ([A-Z0-9]+(?: +[A-Z0-9]+)*?),\s*penalty=\s*([\d.]+)")


@dataclass
class Parameter:
    """One bonded parameter as written in the file.

    values: bond (kb, b0); angle (ktheta, theta0, kub or None, s0 or None); dihedral ((kchi, n, delta), ...) with
    one term per file line; improper (kpsi, n, psi0). annotation: (source types, penalty) from a
    `from <types>, penalty= <p>` comment, i.e. the program's own analogy choice when the parameter was added.
    """
    kind: str
    types: tuple
    values: tuple
    comment: str = ""
    line: int = 0
    annotation: tuple = None


@dataclass
class ParameterFile:
    masses: dict = field(default_factory=dict)            # type -> mass
    bonded: dict = field(default_factory=lambda: {k: {} for k in N_ATOMS})   # kind -> {types: Parameter}, file order
    nonbonded: dict = field(default_factory=dict)         # type -> (eps, rmin/2, eps14 or None, rmin14/2 or None)
    nbfix: dict = field(default_factory=dict)             # (type, type) -> (emin, rmin)


def _annotation(comment, n):
    m = _ANNOTATION.search(comment)
    if not m:
        return None
    src = tuple(m.group(1).split())
    return (src, float(m.group(2))) if len(src) == n else None


def read_parameters(path):
    """Read a CHARMM parameter file (CGenFF layout): ATOMS, BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED, NBFIX."""
    pf = ParameterFile()
    section = None
    with open(path, errors="ignore") as fh:
        for lineno, line in enumerate(fh, 1):
            data, _, comment = line.partition("!")
            words = data.split()
            if not words:
                continue
            head = words[0].upper()
            if head in ("ATOMS", *BONDED, "NONBONDED", "NBFIX", "HBOND", "CMAP", "END"):
                section = head
                continue
            if section == "ATOMS" and head == "MASS":
                pf.masses[words[2]] = float(words[3])
            elif section in BONDED:
                kind = BONDED[section]
                n = N_ATOMS[kind]
                types, nums = tuple(words[:n]), words[n:]
                if kind == "bond":
                    values = (float(nums[0]), float(nums[1]))
                elif kind == "angle":
                    values = (float(nums[0]), float(nums[1]),
                              float(nums[2]) if len(nums) > 2 else None, float(nums[3]) if len(nums) > 3 else None)
                elif kind == "dihedral":
                    values = ((float(nums[0]), int(nums[1]), float(nums[2])),)
                else:
                    values = (float(nums[0]), int(nums[1]), float(nums[2]))
                table = pf.bonded[kind]
                if kind == "dihedral" and types in table:     # further terms of a multi-term dihedral
                    table[types].values += values
                    continue
                if types in table:
                    continue
                table[types] = Parameter(kind, types, values, comment.strip(), lineno, _annotation(comment, n))
            elif section == "NONBONDED" and len(words) >= 4 and words[0] not in ("cutnb",):
                try:
                    nums = [float(x) for x in words[1:7]]
                except ValueError:
                    continue
                pf.nonbonded[words[0]] = (nums[1], nums[2], nums[4] if len(nums) > 4 else None,
                                          nums[5] if len(nums) > 5 else None)
            elif section == "NBFIX" and len(words) >= 4:
                pf.nbfix[(words[0], words[1])] = (float(words[2]), float(words[3]))
    return pf

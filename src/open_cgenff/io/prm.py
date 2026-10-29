"""Reader for the charge-increment field of CGenFF 5.0 parameter files (par_all36_cgenff_increments.prm).

Each bond/angle/dihedral line may carry `!<% d1 [d2 [d3]] %>!`: one increment per bond of the parameter, in the
parameter's atom order; `N/A` = no increment (the term does not contribute to charges). Multi-term dihedrals repeat
the same increments on every line; the first line counts.
"""
import re

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

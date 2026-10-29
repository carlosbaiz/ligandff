"""Parameter database (SPEC 5.8): bonded parameters of par_all36_cgenff.prm with exact lookup and canonical order.

Canonical ordering, from the parameter-file header: element priority C < N < O < P < S < halogens (by Z) < other
elements (by Z) < H, alphabetical within an element. Bonds are written lowest priority first; angles have 1 and 3
swapped if 3 < 1; dihedrals are reversed if 3 < 2 (if 2 = 3, if 4 < 1); impropers are not reordered (first atom
central). Within a section, lines are sorted on the atoms (bonds 1,2; angles 2,1,3; dihedrals 2,3,1,4; impropers
1,4,2,3). The analogy search scans candidates in this collated order, so ties go to the first in that order.
"""
from ..io.prm import N_ATOMS, read_parameters
from ..paths import PARAMETER_FILE, toppar_dir

KINDS = tuple(N_ATOMS)
_SECTION_KEY = {"bond": (0, 1), "angle": (1, 0, 2), "dihedral": (1, 2, 0, 3), "improper": (0, 3, 1, 2)}
# (atomic number, element) for the masses in the parameter file; anything else (lone pairs, unknown) sorts as "misc"
_ELEMENTS = {1.008: (1, "H"), 10.811: (5, "B"), 12.011: (6, "C"), 14.007: (7, "N"), 15.999: (8, "O"),
             18.998: (9, "F"), 22.990: (11, "Na"), 26.982: (13, "Al"), 30.974: (15, "P"), 32.06: (16, "S"),
             35.453: (17, "Cl"), 40.08: (20, "Ca"), 78.96: (34, "Se"), 79.904: (35, "Br"), 126.904: (53, "I")}
_HALOGENS = {"F", "Cl", "Br", "I"}
_GROUP = {"C": 0, "N": 1, "O": 2, "P": 3, "S": 4}


def element_of_mass(mass):
    for m, el in _ELEMENTS.items():
        if abs(mass - m) < 0.01:
            return el
    return (999, "X")


def element_priority(z, symbol):
    """Sort key of an element: C < N < O < P < S < halogens (by Z) < other elements (by Z) < H."""
    if symbol in _GROUP:
        return (_GROUP[symbol], 0)
    if symbol in _HALOGENS:
        return (5, z)
    if symbol == "H":
        return (7, 0)
    return (6, z)


class ParameterDatabase:
    """Bonded parameters by kind ('bond', 'angle', 'dihedral', 'improper') plus nonbonded parameters."""

    def __init__(self, pfile):
        self.file = pfile
        self.elements = {t: element_of_mass(m) for t, m in pfile.masses.items()}
        self._priority = {t: (element_priority(*el), t) for t, el in self.elements.items()}
        self.params = {kind: dict(pfile.bonded[kind]) for kind in KINDS}
        self._collated = {}

    @classmethod
    def from_prm(cls, path):
        return cls(read_parameters(path))

    @classmethod
    def from_toppar(cls, path=None):
        return cls.from_prm(toppar_dir(path) / PARAMETER_FILE)

    @property
    def nonbonded(self):
        return self.file.nonbonded

    def priority(self, atom_type):
        return self._priority.get(atom_type, ((6, 999), atom_type))

    def canonical(self, kind, types):
        """Types in the parameter file's canonical orientation."""
        types = tuple(types)
        p = self.priority
        if kind == "bond":
            return types if p(types[0]) <= p(types[1]) else types[::-1]
        if kind == "angle":
            return types if p(types[0]) <= p(types[2]) else types[::-1]
        if kind == "dihedral":
            a, b, c, d = map(p, types)
            if c < b or (c == b and d < a):
                return types[::-1]
            return types
        return types

    def sort_key(self, kind, types):
        return tuple(self.priority(types[i]) for i in _SECTION_KEY[kind])

    def collated(self, kind):
        """Parameters of one kind in collated order (SPEC 5.8); the analogy search scans them in this order."""
        if kind not in self._collated:
            self._collated[kind] = sorted(self.params[kind].values(), key=lambda prm: self.sort_key(kind, prm.types))
        return self._collated[kind]

    def lookup(self, kind, types):
        """(Parameter, orientation) with orientation +1 if found as given, -1 if reversed; (None, 0) if missing."""
        types = tuple(types)
        table = self.params[kind]
        if types in table:
            return table[types], 1
        if types[::-1] in table:
            return table[types[::-1]], -1
        return None, 0

    def __contains__(self, item):
        kind, types = item
        return self.lookup(kind, types)[0] is not None

    def without(self, kind, types):
        """A copy without one parameter (leave-one-out tests and benchmarks)."""
        db = ParameterDatabase.__new__(ParameterDatabase)
        db.file, db.elements, db._priority = self.file, self.elements, self._priority
        db.params = {k: dict(v) for k, v in self.params.items()}
        prm, _ = self.lookup(kind, types)
        if prm is not None:
            del db.params[kind][prm.types]
        db._collated = {k: v for k, v in self._collated.items() if k != kind}
        return db

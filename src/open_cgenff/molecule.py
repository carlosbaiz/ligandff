"""Molecule containers: InputMolecule (reader output, SPEC 5.1) and TypedMolecule (typing -> parameters/charges)."""
from dataclasses import dataclass, field

AROMATIC = "ar"                                   # input bond order of unknown resonance (mol2 'ar', RDKit aromatic)


@dataclass
class InputMolecule:
    """Molecule as read, in input atom order (resonance and altnum can depend on it, SPEC 5.1).

    bonds: (i, j, order) with order 1, 2, 3 or AROMATIC. formal_charges: per atom, or None when the input format
    does not carry them (mol2); perception then derives them from valences.
    """
    name: str
    atom_names: list
    elements: list
    bonds: list
    formal_charges: list = None
    coords: list = None

    def __post_init__(self):
        n = len(self.elements)
        if len(self.atom_names) != n or (self.formal_charges is not None and len(self.formal_charges) != n):
            raise ValueError("atom_names, elements and formal_charges must have one entry per atom")
        seen = set()
        for i, j, order in self.bonds:
            if not (0 <= i < n and 0 <= j < n) or i == j or frozenset((i, j)) in seen:
                raise ValueError(f"bad or duplicate bond ({i}, {j})")
            if order not in (1, 2, 3, AROMATIC):
                raise ValueError(f"bad bond order {order!r} for ({i}, {j})")
            seen.add(frozenset((i, j)))

    @property
    def n_atoms(self):
        return len(self.elements)

    @classmethod
    def from_rdkit(cls, mol, name=""):
        """Aromatic RDKit bonds become AROMATIC (resolved by perception); formal charges are taken from RDKit."""
        from rdkit import Chem
        orders = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
        bonds = []
        for b in mol.GetBonds():
            order = AROMATIC if b.GetIsAromatic() else orders.get(b.GetBondType())
            if order is None:
                raise ValueError(f"unsupported bond type {b.GetBondType()}")
            bonds.append((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), order))
        atoms = list(mol.GetAtoms())
        return cls(name, [f"{a.GetSymbol()}{a.GetIdx() + 1}" for a in atoms], [a.GetSymbol() for a in atoms], bonds,
                   [a.GetFormalCharge() for a in atoms])


@dataclass
class TypedMolecule:
    """Atom types and typing-rule formal charges of one molecule.

    formal_charges are the integer charges placed by typing-rule `charge` actions (SPEC 5.7: e.g. -1 on the
    carboxylate carbon, +1 on an ammonium N), not the resonance formal charges of perception (SPEC 5.3).
    bonds are 0-based atom index pairs; impropers are 4-tuples from typing-rule `impr` actions.
    """
    names: list
    types: list
    bonds: list
    formal_charges: list = None
    impropers: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def __post_init__(self):
        if self.formal_charges is None:
            self.formal_charges = [0] * len(self.types)
        if not len(self.names) == len(self.types) == len(self.formal_charges):
            raise ValueError("names, types and formal_charges must have one entry per atom")
        n = len(self.types)
        for i, j in self.bonds:
            if not (0 <= i < n and 0 <= j < n) or i == j:
                raise ValueError(f"bad bond ({i}, {j}) for {n} atoms")

    @property
    def net_charge(self):
        return sum(self.formal_charges)

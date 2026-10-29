"""Typed molecule: the hand-off from typing (SPEC 5.5-5.7) to parameter and charge assignment."""
from dataclasses import dataclass, field


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

"""Perception (SPEC 5.2-5.4): rings -> resonance -> ring classes, on an InputMolecule.

    p = perceive(InputMolecule.from_rdkit(mol))          # or io.mol2.read_mol2(path)[0]
    p.order(i, j), p.charges[i], p.atom_rings[i], p.ring_class[ring], p.in_ring(i, j), p.warnings
"""
from dataclasses import dataclass, field

from ..molecule import InputMolecule
from .aromaticity import classify_rings
from .resonance import COEFFICIENTS, COEFFICIENTS_PRE_097, ResonanceResult, resolve
from .rings import atom_rings, canonical_ring, find_rings, ring_bonds


@dataclass
class Perception:
    molecule: InputMolecule
    neighbors: list                   # per atom, ascending atom index
    orders: dict                      # frozenset((i, j)) -> resolved bond order 1/2/3
    charges: list                     # resonance formal charges (not the typing-rule charges of SPEC 5.7)
    rings: list                       # canonical tuples of all rings with < 8 atoms, sorted by (size, atoms)
    ring_class: dict                  # ring -> 'arom' | 'sp2' | 'sp3' | 'mixed'
    atom_rings: list                  # per atom: [smallest ring], or the three smallest for bridgeheads
    bridgehead: list
    resonance: ResonanceResult
    warnings: list = field(default_factory=list)

    def __post_init__(self):
        self._ring_bonds = set().union(*map(ring_bonds, self.rings)) if self.rings else set()

    @property
    def n_atoms(self):
        return len(self.neighbors)

    def order(self, i, j):
        return self.orders.get(frozenset((i, j)), 0)

    def valence(self, i):
        """Number of electrons in bonds (P-I 'nb'): sum of bond orders, hydrogens included."""
        return sum(self.order(i, j) for j in self.neighbors[i])

    def in_ring(self, i, j):
        """Bond i-j is endocyclic (in a ring with < 8 atoms)."""
        return frozenset((i, j)) in self._ring_bonds


def perceive(mol, coefficients=COEFFICIENTS):
    n = mol.n_atoms
    neighbors = [[] for _ in range(n)]
    for i, j, _ in mol.bonds:
        neighbors[i].append(j); neighbors[j].append(i)
    neighbors = [sorted(x) for x in neighbors]
    rings = find_rings(neighbors)
    res = resolve(n, mol.elements, neighbors, mol.bonds, mol.formal_charges, rings, coefficients)
    classes, _, converged = classify_rings(rings, neighbors, mol.elements,
                                           lambda a, b: res.orders.get(frozenset((a, b)), 0), res.charges)
    per_atom, bridgehead = atom_rings(n, neighbors, rings)
    warnings = list(res.warnings)
    if not converged:
        warnings.append("aromaticity perception did not converge")
    return Perception(mol, neighbors, res.orders, res.charges, rings, classes, per_atom, bridgehead, res, warnings)


__all__ = ["Perception", "perceive", "find_rings", "canonical_ring", "classify_rings", "resolve", "COEFFICIENTS",
           "COEFFICIENTS_PRE_097"]

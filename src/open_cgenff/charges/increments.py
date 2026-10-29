"""Charge assignment by bond-charge increments (P-II; SPEC 5.11). Reproduces program 5.0 output exactly on B0.

q_i = q0_i - sum_j beta_ij, extended so that each bond parameter carries one increment, each angle two and each
dihedral three. For an instance A-B(-C(-D)) in parameter-file orientation with increments (d1[, d2[, d3]]), increment
k moves charge along bond k: q_k -= d_k, q_(k+1) += d_k. An instance that matches the parameter in reverse order is
reversed first. Impropers carry no increments; `N/A` increments contribute nothing.

Increments for parameters missing from CGenFF come from an analog parameter supplied by the caller (P-II: found with
the nonbonded penalty matrix at all positions, so it may differ from the bonded-parameter analog). A dihedral analog
with penalty above 50 contributes no charge and is counted with penalty 50 (P-II, 2012; to confirm on B5).
"""
import collections
import itertools
from dataclasses import dataclass, field

from ..io.prm import read_increments
from ..paths import INCREMENTS_FILE, toppar_dir
from .penalties import atom_penalty

MAX_DIHEDRAL_PENALTY = 50.0
_MISSING = object()


@dataclass(frozen=True)
class Analog:
    """Parameter used in place of a missing one: its types (as in the parameter file) and analogy penalty."""
    source: tuple
    penalty: float


class IncrementTable:
    """Charge increments of all bond, angle and dihedral parameters, keyed by atom types."""

    def __init__(self, table):
        self._table = table

    @classmethod
    def from_prm(cls, path):
        return cls(read_increments(path))

    @classmethod
    def from_toppar(cls, path=None):
        return cls.from_prm(toppar_dir(path) / INCREMENTS_FILE)

    def __contains__(self, types):
        return self.lookup(types)[1] != 0

    def lookup(self, types):
        """(increments or None for N/A, orientation): +1 = as written, -1 = reversed, 0 = no such parameter."""
        types = tuple(types)
        t = self._table[len(types)]
        if types in t:
            return t[types], 1
        rev = types[::-1]
        if rev in t:
            return t[rev], -1
        return _MISSING, 0


@dataclass
class ChargeResult:
    charges: list
    penalties: list                                       # per-atom charge penalties (P-II eq 4)
    missing: list                                         # parameter types with neither increments nor an analog
    analogs_used: dict = field(default_factory=dict)      # missing types (as seen in the molecule) -> Analog

    @property
    def charge_penalty(self):
        """RESI header value: maximum per-atom penalty."""
        return max(self.penalties, default=0.0)


def parameter_instances(n_atoms, bonds):
    """All bond, angle and dihedral instances as atom-index tuples, in a fixed order."""
    adj = [set() for _ in range(n_atoms)]
    for i, j in bonds:
        adj[i].add(j); adj[j].add(i)
    out = [tuple(b) for b in bonds]
    for b in range(n_atoms):
        out += [(a, b, c) for a, c in itertools.combinations(sorted(adj[b]), 2)]
    for b, c in bonds:
        out += [(a, b, c, d) for a in adj[b] - {c} for d in adj[c] - {b, a}]
    return out


def _normalize_analogs(analogs):
    out = {}
    for key, val in (analogs or {}).items():
        if not isinstance(val, Analog):
            val = Analog(tuple(val[0]), float(val[1]))
        out[tuple(key)] = val
    return out


def assign_charges(types, bonds, formal_charges, table, analogs=None, max_dihedral_penalty=MAX_DIHEDRAL_PENALTY):
    """Charges and per-atom charge penalties.

    types: atom types; bonds: 0-based index pairs; formal_charges: typing-rule formal charges (SPEC 5.7);
    table: IncrementTable; analogs: {missing types: Analog or (source types, penalty)}, either orientation.
    """
    n = len(types)
    analogs = _normalize_analogs(analogs)
    q = [float(x) for x in formal_charges]
    contrib = collections.defaultdict(list)               # atom -> [(net transfer, penalty)]
    missing, used = set(), {}
    for inst in parameter_instances(n, bonds):
        key = tuple(types[a] for a in inst)
        values, orient = table.lookup(key)
        penalty = 0.0
        if values is _MISSING:
            analog, key_orient = analogs.get(key), 1
            if analog is None:
                analog, key_orient = analogs.get(key[::-1]), -1
            if analog is None:
                missing.add(key); continue
            values, src_orient = table.lookup(analog.source)
            if values is _MISSING:
                missing.add(key); continue
            used[key] = analog
            orient, penalty = key_orient * src_orient, analog.penalty
            if len(inst) == 4 and penalty > max_dihedral_penalty:
                values, penalty = (0.0, 0.0, 0.0), max_dihedral_penalty
        if values is None:
            continue
        atoms = inst if orient == 1 else inst[::-1]
        net = collections.defaultdict(float)
        for k, d in enumerate(values):
            net[atoms[k]] -= d
            net[atoms[k + 1]] += d
        for a, dq in net.items():
            q[a] += dq
            if penalty > 0:
                contrib[a].append((dq, penalty))
    return ChargeResult(q, [atom_penalty(contrib[i]) for i in range(n)], sorted(missing), used)


def charge_molecule(mol, table, analogs=None, **kwargs):
    """assign_charges for a TypedMolecule, using the formal charges from typing."""
    return assign_charges(mol.types, mol.bonds, mol.formal_charges, table, analogs, **kwargs)

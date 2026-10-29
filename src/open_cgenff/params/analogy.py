"""Bonded parameters by analogy (P-II; SPEC 5.9) and parameter assignment for a typed molecule.

For a missing parameter, every existing parameter of the same kind is scored with the total penalty score
TPS = sum over positions of weight x matrix[missing type, candidate type] + bond-group penalties, for the missing
parameter as given and reversed (impropers: also each order of the three outer atoms); the lowest TPS wins, the first
candidate in collated order on ties (SPEC 5.8). Inner positions (bond atoms, angle centre, dihedral 2-3, improper
centre) use the bonded matrix and weight 10, except improper outer atoms (bonded matrix, weight 1); outer positions
(angle and dihedral ends) use the nonbonded matrix. Charge-increment analogs (SPEC 5.11) use the nonbonded matrix at
every position, same weights and bond groups, and only candidates that carry increments.
"""
import itertools
from dataclasses import dataclass, field

from ..charges.increments import Analog, parameter_instances
from .penalties import bond_weights, position_roles

KIND_OF_SIZE = {2: "bond", 3: "angle"}
INNER = {"bond": (0, 1), "angle": (1,), "dihedral": (1, 2), "improper": (0,)}      # inner positions
INNER_BONDS = {"bond": (0,), "angle": (), "dihedral": (1,), "improper": ()}       # virtual bonds within them


@dataclass(frozen=True)
class AnalogResult:
    """Best analog of a missing parameter."""
    kind: str
    missing: tuple            # missing types as passed to the search
    parameter: object         # io.prm.Parameter of the source
    aligned: tuple            # source types aligned position by position with `missing`
    penalty: float

    @property
    def source(self):
        """Source types as written in the parameter file."""
        return self.parameter.types


class AnalogyEngine:
    """Analogy search over a ParameterDatabase with a PenaltyModel; results are cached."""

    def __init__(self, db, model):
        self.db = db
        self.model = model
        self._cache = {}
        self._by_inner = {}

    # ---- scoring ------------------------------------------------------------------------------------------------
    def roles(self, kind, mode):
        roles = position_roles(kind)
        if mode == "charge":
            roles = tuple(("nonbonded", w) for _, w in roles)
        return roles

    def tps(self, kind, missing, candidate, mode="bonded"):
        """TPS of `candidate` for `missing`, both in the given orientation (no reversal)."""
        total = 0.0
        for (mat, w), a, b in zip(self.roles(kind, mode), missing, candidate):
            if a != b:
                total += w * self.model.matrix(mat).cost(a, b)
        groups = self.model.groups
        for k, w in enumerate(bond_weights(kind)):
            total += w * groups.penalty((missing[k], missing[k + 1]), (candidate[k], candidate[k + 1]))
        return total

    @staticmethod
    def orientations(kind, n):
        """Index orders of the missing parameter that are evaluated: as given, reversed; impropers: centre first and
        every order of the outer atoms (as given first)."""
        if kind == "improper":
            return [(0, *p) for p in itertools.permutations(range(1, n))]
        return [tuple(range(n)), tuple(range(n - 1, -1, -1))]

    # ---- search -------------------------------------------------------------------------------------------------
    def _index(self, kind):
        """Candidates grouped by their inner-position types: {inner types: [(collated index, Parameter)]}."""
        if kind not in self._by_inner:
            groups = {}
            for n, prm in enumerate(self.db.collated(kind)):
                groups.setdefault(tuple(prm.types[i] for i in INNER[kind]), []).append((n, prm))
            self._by_inner[kind] = groups
        return self._by_inner[kind]

    def find(self, kind, types, mode="bonded", accept=None):
        """Best AnalogResult for missing `types` (None if no candidate). accept: optional predicate on candidates;
        mode 'bonded' (parameter values) or 'charge' (increments, nonbonded matrix everywhere).

        Equivalent to scoring every candidate in collated order and keeping the first lowest TPS; candidates are
        grouped by inner-position types and groups whose inner cost alone exceeds the best TPS are skipped."""
        key = (kind, tuple(types), mode, accept)
        if key in self._cache:
            return self._cache[key]
        types = tuple(types)
        perms = self.orientations(kind, len(types))
        orients = [tuple(types[i] for i in p) for p in perms]
        roles = self.roles(kind, mode)
        mats = [self.model.matrix(m) for m, _ in roles]
        weights = [w for _, w in roles]
        bweights = bond_weights(kind)
        groups = self.model.groups
        inner, inner_bonds = INNER[kind], INNER_BONDS[kind]
        cost_cache = [[{} for _ in roles] for _ in orients]

        def cost(oi, i, b):
            a = orients[oi][i]
            if a == b:
                return 0.0
            c = cost_cache[oi][i].get(b)
            if c is None:
                c = cost_cache[oi][i][b] = weights[i] * mats[i].cost(a, b)
            return c

        missing_groups = [[groups.membership(o[k], o[k + 1]) for k in range(len(bweights))] for o in orients]

        def group_penalty(oi, k, cand):
            cg = groups.membership(cand[k], cand[k + 1])
            mg = missing_groups[oi][k]
            if cg == mg:
                return 0.0
            return bweights[k] * sum(pen for (pen, _), x, y in zip(groups.groups, mg, cg) if x != y)

        index = self._index(kind)
        bounds = []                                       # (lower bound, orientation, inner types)
        for oi in range(len(orients)):
            for itypes in index:
                full = [None] * len(types)
                for i, t in zip(inner, itypes):
                    full[i] = t
                lb = sum(cost(oi, i, t) for i, t in zip(inner, itypes))
                lb += sum(group_penalty(oi, k, full) for k in inner_bonds)
                bounds.append((lb, oi, itypes))
        bounds.sort(key=lambda x: x[0])
        best = None                                       # (score, collated index, orientation, Parameter)
        for lb, oi, itypes in bounds:
            if best is not None and lb > best[0] + 1e-9:
                break
            for n, prm in index[itypes]:
                if accept is not None and not accept(prm):
                    continue
                cand = prm.types
                score = sum(cost(oi, i, b) for i, b in enumerate(cand))
                score += sum(group_penalty(oi, k, cand) for k in range(len(bweights)))
                if (best is None or score < best[0] - 1e-9
                        or (abs(score - best[0]) <= 1e-9 and (n, oi) < (best[1], best[2]))):
                    best = (score, n, oi, prm)
        result = None
        if best is not None:
            score, _, oi, prm = best
            aligned = [None] * len(types)
            for j, i in enumerate(perms[oi]):
                aligned[i] = prm.types[j]
            result = AnalogResult(kind, types, prm, tuple(aligned), round(score, 6))
        self._cache[key] = result
        return result


# ---- assignment for a molecule ------------------------------------------------------------------------------------

@dataclass
class AssignedParameter:
    kind: str
    atoms: tuple              # atom indices, oriented like `types`
    types: tuple
    parameter: object         # io.prm.Parameter (the exact one, or the analog's source)
    penalty: float = 0.0
    analog: AnalogResult = None

    @property
    def is_analog(self):
        return self.analog is not None


@dataclass
class ParameterAssignment:
    parameters: list = field(default_factory=list)       # [AssignedParameter], one per instance
    missing: list = field(default_factory=list)          # [(kind, types)] with no candidate at all

    @property
    def param_penalty(self):
        """RESI header value: maximum penalty over assigned parameters."""
        return max((p.penalty for p in self.parameters), default=0.0)

    def analogs(self):
        """Distinct analog parameters: {(kind, types): AnalogResult} in first-use order."""
        out = {}
        for p in self.parameters:
            if p.analog is not None:
                out.setdefault((p.kind, p.types), p.analog)
        return out


def molecule_instances(mol):
    """(kind, atoms) for every bond, angle, dihedral and typing-rule improper of a TypedMolecule."""
    out = []
    for inst in parameter_instances(len(mol.types), mol.bonds):
        out.append((KIND_OF_SIZE.get(len(inst), "dihedral"), tuple(inst)))
    out += [("improper", tuple(imp)) for imp in mol.impropers]
    return out


def _orient(db, kind, atoms, types):
    """Atoms and types in the parameter file's canonical order (impropers: as given)."""
    canon = db.canonical(kind, types)
    if canon != tuple(types):
        return tuple(atoms[::-1]), canon
    return tuple(atoms), tuple(types)


def assign_parameters(mol, engine):
    """Exact parameters or analogs for every bonded term of a TypedMolecule (SPEC 5.9)."""
    db = engine.db
    out = ParameterAssignment()
    for kind, atoms in molecule_instances(mol):
        types = tuple(mol.types[a] for a in atoms)
        if kind == "improper":
            prm, order = _find_improper(db, types)
            if prm is not None:
                atoms = tuple(atoms[i] for i in order)
                out.parameters.append(AssignedParameter(kind, atoms, prm.types, prm))
                continue
        else:
            atoms, types = _orient(db, kind, atoms, types)
            prm, orient = db.lookup(kind, types)
            if prm is not None:
                if orient == -1:
                    atoms, types = atoms[::-1], types[::-1]
                out.parameters.append(AssignedParameter(kind, atoms, types, prm))
                continue
        res = engine.find(kind, types)
        if res is None:
            out.missing.append((kind, types))
            continue
        out.parameters.append(AssignedParameter(kind, atoms, types, res.parameter, res.penalty, res))
    return out


def _find_improper(db, types):
    """Exact improper: (Parameter, atom order), centre first, trying every order of the outer atoms (as given first).
    Not reversed: the file's first atom is the central one."""
    table = db.params["improper"]
    for perm in itertools.permutations((1, 2, 3)):
        order = (0, *perm)
        prm = table.get(tuple(types[i] for i in order))
        if prm is not None:
            return prm, order
    return None, None


def charge_analogs(mol, engine, table):
    """{missing types: charges.Analog} for the bonds, angles and dihedrals of `mol` without charge increments in
    `table` (an IncrementTable): analogs found with the nonbonded matrix at all positions (SPEC 5.11)."""
    accept = _increment_filter(table)
    out = {}
    for kind, atoms in molecule_instances(mol):
        if kind == "improper":
            continue
        types = tuple(mol.types[a] for a in atoms)
        if types in table or types in out or types[::-1] in out:
            continue
        res = engine.find(kind, types, mode="charge", accept=accept)
        if res is not None:
            out[types] = Analog(res.aligned, res.penalty)
    return out


_FILTERS = {}


def _increment_filter(table):
    """Candidate predicate: the parameter has numeric increments in `table` (one cached function per table so the
    analogy cache can key on it)."""
    f = _FILTERS.get(id(table))
    if f is None or f[0] is not table:
        def accept(prm, _table=table):
            values, orient = _table.lookup(prm.types)
            return orient != 0 and values is not None
        f = _FILTERS[id(table)] = (table, accept)
    return f[1]

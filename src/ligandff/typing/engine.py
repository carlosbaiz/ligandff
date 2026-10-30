"""Interpreter for the typing rule language (P-I; SPEC 5.5) on a Perception.

Typing of each atom starts at `cat main`; in each category the first rule whose conditions all hold fires. `ne`
matches its groups greedily and in order against the not-yet-used neighbours, lowest atom index first, without
backtracking. Inside a group evaluated on a neighbour, `bo`/`inring` refer to the bond just traversed, and nested
`ne` does not go back over that bond (interpretation: otherwise `ne ( ne ( el C ) )` would trivially match the atom
being typed; `self` then detects ring closures). A ring matched by a ring keyword is consumed for the rest of that
group, so `arom 6 arom 6` needs two aromatic 6-rings. Ring keywords look at the atom's ring list (smallest ring, or
the three smallest for bridgeheads, SPEC 5.2); `inring` and `rings` as noted in rules.py.
"""
from dataclasses import dataclass, field
from fnmatch import fnmatchcase

from .rules import RuleSet

HALOGENS = frozenset(("F", "Cl", "Br", "I"))
MAX_SUB_DEPTH = 64


class TypingError(Exception):
    pass


@dataclass
class AtomTyping:
    type: str
    charge: int = 0                   # typing-rule formal charge (SPEC 5.7)
    impr: bool = False
    altnum: bool = False
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    trace: list = field(default_factory=list)       # rules that fired, main first


@dataclass
class TypingResult:
    types: list
    charges: list
    impropers: list                   # (center, n1, n2, n3)
    atoms: list                       # AtomTyping per atom
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)


class _Evaluator:
    def __init__(self, rules, perception):
        p = perception
        self.rules, self.p = rules, p
        self.el = p.molecule.elements
        self.nbrs = p.neighbors
        self.valence = [p.valence(i) for i in range(p.n_atoms)]
        self.n_h = [sum(self.el[j] == "H" for j in p.neighbors[i]) for i in range(p.n_atoms)]
        self.bond_rings = {}
        for r in p.rings:
            for k in range(len(r)):
                self.bond_rings.setdefault(frozenset((r[k], r[(k + 1) % len(r)])), []).append(r)
        # fusion atom: bridgehead (SPEC 5.2) two of whose rings share at least one more atom (spiro atoms are not)
        # and that has >= 3 neighbours in the two rings together
        self.fusion = [p.bridgehead[i] and any(len(set(r1) & set(r2)) >= 2
                                               and sum(j in set(r1) | set(r2) for j in p.neighbors[i]) >= 3
                                               for x, r1 in enumerate(p.atom_rings[i]) for r2 in p.atom_rings[i][x + 1:])
                       for i in range(p.n_atoms)]
        self.typed = [None] * p.n_atoms
        self.active = set()
        self.root = None

    def _ring_ok(self, r, cls, size):
        return (cls is None or self.p.ring_class[r] == cls) and (size is None or len(r) == size)

    def match(self, group, a, came, consumed):
        return all(self.cond(c, a, came, consumed) for c in group)

    def cond(self, c, a, came, consumed):
        k = c[0]
        if k == "el":
            return self.el[a] == c[1]
        if k == "elha":
            return self.el[a] in HALOGENS
        if k == "elos":
            return self.el[a] in ("O", "S")
        if k == "nb":
            return self.valence[a] == c[1]
        if k == "nh":
            return self.n_h[a] == c[1]
        if k == "deg":
            return len(self.nbrs[a]) == c[1]
        if k == "fc":
            return self.p.charges[a] == c[1]
        if k == "rings":
            return len(self.p.atom_rings[a]) == c[1]
        if k == "ringkw":
            for r in self.p.atom_rings[a]:
                if r not in consumed and self._ring_ok(r, c[1], c[2]):
                    consumed.add(r)
                    return True
            return False
        if k == "self":
            return a == self.root
        if k == "fusion":
            return self.fusion[a]
        if k == "ne":
            used = set()
            for g in c[1]:
                for b in self.nbrs[a]:
                    if b != came and b not in used and self.match(g, b, a, set()):
                        used.add(b)
                        break
                else:
                    return False
            return True
        if k == "bo":
            return came is not None and self.p.order(came, a) == c[1]
        if k == "inring":
            if came is None:
                return False
            return any(self._ring_ok(r, c[1], c[2]) for r in self.bond_rings.get(frozenset((came, a)), ()))
        if k == "not":
            return not self.match(c[1], a, came, set(consumed))
        if k == "or":
            for g in c[1]:
                trial = set(consumed)
                if self.match(g, a, came, trial):
                    consumed |= trial
                    return True
            return False
        if k == "type":
            return fnmatchcase(self.type_atom(a).type, c[1])
        if k == "rmember":
            return any(self.match(c[3], m, None, set())
                       for r in self.p.atom_rings[a] if self._ring_ok(r, c[1], c[2]) for m in r)
        raise TypingError(f"unknown condition {k}")

    def type_atom(self, i):
        if self.typed[i] is not None:
            return self.typed[i]
        if i in self.active:
            raise TypingError(f"atom {i + 1}: circular 'type' dependency")
        self.active.add(i)
        saved_root, self.root = self.root, i
        try:
            result = self._run(i)
        finally:
            self.root = saved_root
            self.active.discard(i)
        self.typed[i] = result
        return result

    def _run(self, i):
        out = AtomTyping(type=None)
        cat = "main"
        for _ in range(MAX_SUB_DEPTH):
            for rule in self.rules.categories[cat]:
                if self.match(rule.conditions, i, None, set()):
                    break
            else:
                out.errors.append(f"no rule matches in category {cat}")
                out.type = "UNK"
                return out
            out.trace.append(rule)
            if rule.charge is not None:
                out.charge = rule.charge
            out.impr |= rule.impr
            out.altnum |= rule.altnum
            out.warnings += rule.warn
            out.errors += rule.err
            if rule.action == "typ":
                out.type = rule.target
                return out
            cat = rule.target
        raise TypingError(f"atom {i + 1}: more than {MAX_SUB_DEPTH} nested 'sub' jumps")


def _resolve_altnum(types, atoms, perception):
    """Replace '?' by 1/2: atoms joined by a double bond share the digit, by a single bond they alternate (SPEC 5.5).
    Each conjugated chain is numbered from its lowest-index atom, which gets 1."""
    todo = {i for i, a in enumerate(atoms) if a.altnum and "?" in a.type}
    digit, warnings = {}, []
    for start in sorted(todo):
        if start in digit:
            continue
        digit[start] = 1
        stack = [start]
        while stack:
            i = stack.pop()
            for j in perception.neighbors[i]:
                if j not in todo:
                    continue
                want = digit[i] if perception.order(i, j) >= 2 else 3 - digit[i]
                if j not in digit:
                    digit[j] = want
                    stack.append(j)
                elif digit[j] != want:
                    warnings.append(f"altnum: inconsistent numbering at atoms {i + 1}-{j + 1}")
    out = list(types)
    for i, d in digit.items():
        out[i] = out[i].replace("?", str(d), 1)
    return out, warnings


def type_perception(rules, perception):
    """Type every atom of a perceived molecule."""
    if not isinstance(rules, RuleSet):
        raise TypeError("rules must be a RuleSet")
    ev = _Evaluator(rules, perception)
    atoms = [ev.type_atom(i) for i in range(perception.n_atoms)]
    types, warnings = _resolve_altnum([a.type for a in atoms], atoms, perception)
    impropers = [(i, *perception.neighbors[i]) for i, a in enumerate(atoms) if a.impr and len(perception.neighbors[i]) == 3]
    for i, a in enumerate(atoms):
        warnings += [f"atom {i + 1} ({types[i]}): {w}" for w in a.warnings]
    errors = [f"atom {i + 1} ({types[i]}): {e}" for i, a in enumerate(atoms) for e in a.errors]
    return TypingResult(types, [a.charge for a in atoms], impropers, atoms, warnings, errors)

"""Penalty rules and the substitution matrices for the analogy search (SPEC 5.9-5.10).

Two asymmetric matrices give the cost of substituting the atom type of a missing parameter (row) by the type of a
candidate parameter (column): the bonded matrix for inner positions, the nonbonded matrix for outer positions.

Provisional matrix (M4), looked up in this order:
1. identical types: 0;
2. the hierarchical `pri/alt/up` categories of the penalty rule file (P-II Table 1, bonded matrix only);
3. elements observed in the parameter file: an annotated parameter `from <source>, penalty= <tps>` that differs from
   its source at exactly one position, with no bond-group difference, gives one element = tps / position weight
   (FINDINGS §4); the most frequent value wins (ties: the lowest);
4. a fallback by type similarity: the median of the observed elements of the same matrix with the same relation
   (`_relation`: different elements, same element, same family).

Fitted matrix (M5, `fitted_model`, the default): the provisional matrix learned from B4-train annotations only, with
altnum twins at 0, elements solved from multi-substitution annotations, and the fallback times FALLBACK_SCALE
(unobserved substitutions are expensive; FINDINGS §10).
"""
import collections
import re
import statistics
from dataclasses import dataclass, field
from pathlib import Path

RULES_DIR = Path(__file__).resolve().parent / "rules"
PENALTY_RULES_FILE = RULES_DIR / "penalties_v5.rules"
MATRICES = ("bonded", "nonbonded")
INNER_WEIGHT = 10.0


class PenaltyRuleError(ValueError):
    pass


# ---- rule file -------------------------------------------------------------------------------------------------------

@dataclass
class Entry:
    name: str
    is_type: bool
    pri: float = 0.0
    alts: dict = field(default_factory=dict)
    up: float = None


@dataclass
class Category:
    name: str
    matrix: str
    entries: dict = field(default_factory=dict)          # name -> Entry, file order


@dataclass
class BondGroup:
    penalty: float
    types: frozenset


@dataclass
class PenaltyRules:
    categories: dict = field(default_factory=dict)       # name -> Category
    bond_groups: list = field(default_factory=list)      # [BondGroup], file order

    @classmethod
    def from_file(cls, path=PENALTY_RULES_FILE):
        return parse_penalty_rules(Path(path).read_text(), str(path))

    # cost from the hierarchy -------------------------------------------------------------------------------------
    def _parents(self, matrix):
        """entry name -> (category containing it) for one matrix."""
        out = {}
        for cat in self.categories.values():
            if cat.matrix == matrix:
                for name in cat.entries:
                    out[name] = cat
        return out

    def hierarchy_cost(self, matrix, a, b):
        """Substitution cost a -> b from the categories, or None if the categories do not connect a and b."""
        if a == b:
            return 0.0
        parent = self._parents(matrix)
        if a not in parent or b not in parent:
            return None

        def path(x):                                      # [(entry, category)] from the type upwards
            out = []
            while x in parent:
                cat = parent[x]
                out.append((cat.entries[x], cat))
                x = cat.name
            return out

        pa, pb = path(a), path(b)
        cats_b = {cat.name: i for i, (_, cat) in enumerate(pb)}
        climb = 0.0
        for entry, cat in pa:
            if cat.name in cats_b:
                target, _ = pb[cats_b[cat.name]]
                alt = entry.alts.get(target.name)
                if alt is None:
                    return None
                descend = sum(e.pri for e, _ in pb[:cats_b[cat.name]])
                return climb + alt + descend
            if entry.up is None:
                return None
            climb += entry.up
        return None

    def types(self, matrix):
        return {name for cat in self.categories.values() if cat.matrix == matrix
                for name, e in cat.entries.items() if e.is_type}


_NUMBER = re.compile(r"^-?\d+(\.\d*)?$")


def _number(word, where):
    if not _NUMBER.match(word):
        raise PenaltyRuleError(f"{where}: expected a number, got {word!r}")
    return float(word)


def parse_penalty_rules(text, source="<string>"):
    rules = PenaltyRules()
    statements = []                                       # (line number, words); indented lines continue `bgrp`
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#")[0]
        words = line.split()
        if not words:
            continue
        if raw[:1].isspace() and statements and statements[-1][1][0] == "bgrp":
            statements[-1][1].extend(words)
        else:
            statements.append((lineno, words))
    matrix, cat, subs = "bonded", None, []
    for lineno, words in statements:
        where = f"{source}:{lineno}"
        head = words[0]
        if head == "matrix":
            if len(words) != 2 or words[1] not in MATRICES:
                raise PenaltyRuleError(f"{where}: matrix bonded|nonbonded")
            matrix = words[1]
        elif head == "cat":
            if len(words) != 2:
                raise PenaltyRuleError(f"{where}: cat NAME")
            if words[1] in rules.categories:
                raise PenaltyRuleError(f"{where}: duplicate category {words[1]}")
            cat = rules.categories[words[1]] = Category(words[1], matrix)
        elif head == "end":
            cat = None
        elif head == "bgrp":
            if len(words) < 3:
                raise PenaltyRuleError(f"{where}: bgrp PENALTY TYPE...")
            rules.bond_groups.append(BondGroup(_number(words[1], where), frozenset(words[2:])))
        elif head in ("typ", "sub"):
            if cat is None:
                raise PenaltyRuleError(f"{where}: {head} outside a category")
            if len(words) < 3 or words[2] != ":":
                raise PenaltyRuleError(f"{where}: {head} NAME : ...")
            entry = Entry(words[1], head == "typ")
            rest = words[3:]
            i = 0
            while i < len(rest):
                key = rest[i]
                if key == "pri" and i + 1 < len(rest):
                    entry.pri = _number(rest[i + 1], where); i += 2
                elif key == "up" and i + 1 < len(rest):
                    entry.up = _number(rest[i + 1], where); i += 2
                elif key == "alt" and i + 2 < len(rest):
                    entry.alts[rest[i + 1]] = _number(rest[i + 2], where); i += 3
                else:
                    raise PenaltyRuleError(f"{where}: unexpected {key!r}")
            if entry.name in cat.entries:
                raise PenaltyRuleError(f"{where}: duplicate entry {entry.name} in {cat.name}")
            cat.entries[entry.name] = entry
            if not entry.is_type:
                subs.append((where, entry.name))
        else:
            raise PenaltyRuleError(f"{where}: unknown statement {head!r}")
    for where, name in subs:
        if name not in rules.categories:
            raise PenaltyRuleError(f"{where}: sub to undefined category {name}")
    for cat in rules.categories.values():
        for entry in cat.entries.values():
            for alt in entry.alts:
                if alt not in cat.entries:
                    raise PenaltyRuleError(f"{source}: {cat.name}/{entry.name}: alt {alt} is not in the category")
    return rules


# ---- bond groups -----------------------------------------------------------------------------------------------------

class BondGroups:
    """Bond-group penalties (P-II, program >= 0.9.6). The first two groups count as one (equivalent)."""

    def __init__(self, groups):
        groups = list(groups)
        merged = []
        if len(groups) >= 2:
            merged.append((groups[0].penalty, (groups[0].types, groups[1].types)))
            groups = groups[2:]
        merged += [(g.penalty, (g.types,)) for g in groups]
        self.groups = merged

    def membership(self, a, b):
        """Tuple of booleans: is the bond a-b in each group (both types in one of the group's type sets)."""
        return tuple(any(a in s and b in s for s in sets) for _, sets in self.groups)

    def penalty(self, missing_bond, candidate_bond):
        """Unweighted group penalty for one virtual bond."""
        if not self.groups:
            return 0.0
        m, c = self.membership(*missing_bond), self.membership(*candidate_bond)
        return sum(pen for (pen, _), x, y in zip(self.groups, m, c) if x != y)


# ---- positions and weights -----------------------------------------------------------------------------------------

def position_roles(kind):
    """Per position: (matrix, weight) for the bonded-parameter search (SPEC 5.9)."""
    b, n = ("bonded", INNER_WEIGHT), ("nonbonded", 1.0)
    return {"bond": (b, b), "angle": (n, b, n), "dihedral": (n, b, b, n),
            "improper": (b, ("bonded", 1.0), ("bonded", 1.0), ("bonded", 1.0))}[kind]


def bond_weights(kind):
    """Weights of the virtual bonds (1/2/3 for bond/angle/dihedral; none for impropers)."""
    return {"bond": (INNER_WEIGHT,), "angle": (INNER_WEIGHT, INNER_WEIGHT),
            "dihedral": (1.0, INNER_WEIGHT, 1.0), "improper": ()}[kind]


# ---- matrices --------------------------------------------------------------------------------------------------------

_FAMILY = re.compile(r"^([A-Z]+?G.?)(.*)$")


def _relation(a, b, elements):
    """('diff',): different elements; ('elem', is_H): same element, different family (element + G + hybridization or
    class letter, e.g. CG3, CG2, HGA, NG2); ('fam', n): same family and n further common name characters (max 2).
    H is kept apart because H types substitute cheaply at bonded positions (P-II: fitted from bond lengths)."""
    ea, eb = elements.get(a), elements.get(b)
    if ea is None or eb is None or ea != eb:
        return ("diff",)
    ma, mb = _FAMILY.match(a), _FAMILY.match(b)
    if not ma or not mb or ma.group(1) != mb.group(1):
        return ("elem", ea == "H")
    ra, rb = ma.group(2), mb.group(2)
    n = 0
    while n < min(len(ra), len(rb), 2) and ra[n] == rb[n]:
        n += 1
    return ("fam", n)


RELATIONS = (("fam", 2), ("fam", 1), ("fam", 0), ("elem", True), ("elem", False), ("diff",))   # closest first


class SubstitutionMatrix:
    def __init__(self, name, explicit, fallback, rules=None, elements=None):
        self.name = name
        self.explicit = explicit                          # (a, b) -> observed element
        self.fallback = fallback                          # relation -> value
        self.rules = rules
        self.elements = elements or {}
        self._cache = {}

    def source(self, a, b):
        """Where cost(a, b) comes from: 'identity', 'rules', 'observed' or 'fallback'."""
        if a == b:
            return "identity"
        if self.rules is not None and self.rules.hierarchy_cost(self.name, a, b) is not None:
            return "rules"
        return "observed" if (a, b) in self.explicit else "fallback"

    def cost(self, a, b):
        key = (a, b)
        if key not in self._cache:
            src = self.source(a, b)
            if src == "identity":
                v = 0.0
            elif src == "rules":
                v = self.rules.hierarchy_cost(self.name, a, b)
            elif src == "observed":
                v = self.explicit[key]
            else:
                rel = _relation(a, b, self.elements)
                v = self.fallback[rel]
            self._cache[key] = v
        return self._cache[key]


@dataclass
class PenaltyModel:
    bonded: SubstitutionMatrix
    nonbonded: SubstitutionMatrix
    groups: BondGroups

    def matrix(self, name):
        return self.bonded if name == "bonded" else self.nonbonded


def observed_elements(db, groups, use=None):
    """{matrix: {(missing type, source type): [element values]}} from single-substitution annotations.
    use: optional predicate (kind, Parameter) -> bool selecting the annotations to learn from (e.g. B4-train)."""
    obs = {m: collections.defaultdict(list) for m in MATRICES}
    for kind, table in db.params.items():
        roles = position_roles(kind)
        for prm in table.values():
            if prm.annotation is None or (use is not None and not use(kind, prm)):
                continue
            src, tps = prm.annotation
            best = None
            for cand in (src, src[::-1]):
                diffs = [i for i, (x, y) in enumerate(zip(prm.types, cand)) if x != y]
                if best is None or len(diffs) < len(best[0]):
                    best = (diffs, cand)
            diffs, cand = best
            if len(diffs) != 1:
                continue
            if any(groups.penalty((prm.types[k], prm.types[k + 1]), (cand[k], cand[k + 1]))
                   for k in range(len(bond_weights(kind)))):
                continue
            i = diffs[0]
            mat, weight = roles[i]
            obs[mat][(prm.types[i], cand[i])].append(round(tps / weight, 4))
    return obs


def _mode(values):
    counts = collections.Counter(values)
    best = max(counts.values())
    return min(v for v, c in counts.items() if c == best)


def provisional_model(db, rules=None, use=None):
    """The M4 provisional penalty model (module docstring) for a ParameterDatabase; `use` as in observed_elements."""
    rules = rules if rules is not None else PenaltyRules.from_file()
    groups = BondGroups(rules.bond_groups)
    elements = {t: el[1] for t, el in db.elements.items()}
    obs = observed_elements(db, groups, use)
    mats = {}
    for m in MATRICES:
        explicit = {k: _mode(v) for k, v in obs[m].items()}
        by_rel = collections.defaultdict(list)
        for (a, b), v in explicit.items():
            by_rel[_relation(a, b, elements)].append(v)
        median = {rel: statistics.median(v) for rel, v in by_rel.items()}
        median.setdefault(("diff",), max(explicit.values(), default=0.0))
        # unobserved relation: the next farther observed one; within a family, fewer common characters never cost less
        fallback, farther = {}, median[("diff",)]
        for rel in reversed(RELATIONS):
            farther = fallback[rel] = median.get(rel, farther)
        for n in (1, 0):
            fallback[("fam", n)] = max(fallback[("fam", n)], fallback[("fam", n + 1)])
        mats[m] = SubstitutionMatrix(m, explicit, fallback, rules, elements)
    return PenaltyModel(mats["bonded"], mats["nonbonded"], groups)


# ---- fitted model (M5) ---------------------------------------------------------------------------------------------

# altnum twins (SPEC 5.5, 7.3): the 1/2 digit is arbitrary and parameter-symmetric, so substituting one for the other
# costs nothing
ALTNUM_TWINS = (("CG2DC1", "CG2DC2"), ("CG2D1O", "CG2D2O"), ("CG25C1", "CG25C2"), ("CG251O", "CG252O"))
# unobserved substitutions: relation median x this factor (B4-train grid 5/10/20/50: 91.8/92.4/93.0/93.2 % same source)
FALLBACK_SCALE = 20.0
SPLIT_SEED, TEST_FRACTION = 42, 0.2


def annotation_split(db, seed=SPLIT_SEED, test_fraction=TEST_FRACTION):
    """B4 split of the annotated parameters: {"train": [...], "test": [...]} of "kind:T1 T2 ..." names (sorted)."""
    import random
    names = sorted(f"{kind}:{' '.join(p.types)}" for kind, table in db.params.items()
                   for p in table.values() if p.annotation is not None)
    random.Random(seed).shuffle(names)
    n_test = round(len(names) * test_fraction)
    return {"train": sorted(names[n_test:]), "test": sorted(names[:n_test])}


def _source_path(kind, types, source):
    """The source orientation with the fewest substituted positions (as given first)."""
    best = None
    for cand in (tuple(source), tuple(source)[::-1]):
        diffs = [i for i, (a, b) in enumerate(zip(types, cand)) if a != b]
        if best is None or len(diffs) < len(best[0]):
            best = (diffs, cand)
    return best


def solved_elements(db, model, use=None, rounds=10):
    """Elements fixed by annotations with several substitutions once all but one of their elements are known
    (repeated until nothing changes): {(matrix, missing type, source type): value}."""
    equations = []
    for kind, table in db.params.items():
        roles = position_roles(kind)
        for prm in table.values():
            if prm.annotation is None or (use is not None and not use(kind, prm)):
                continue
            src, tps = prm.annotation
            diffs, cand = _source_path(kind, prm.types, src)
            groups = sum(w * model.groups.penalty((prm.types[k], prm.types[k + 1]), (cand[k], cand[k + 1]))
                         for k, w in enumerate(bond_weights(kind)))
            equations.append(([((roles[i][0], prm.types[i], cand[i]), roles[i][1]) for i in diffs], tps - groups))

    def known(e):
        mat, a, b = e
        return e in solved or model.matrix(mat).source(a, b) != "fallback"

    def value(e):
        return solved[e] if e in solved else model.matrix(e[0]).cost(e[1], e[2])

    solved = {}
    for _ in range(rounds):
        new = collections.defaultdict(list)
        for terms, rhs in equations:
            unknown = [(e, w) for e, w in terms if not known(e)]
            if len(unknown) != 1:
                continue
            (e, w), = unknown
            v = round((rhs - sum(w2 * value(e2) for e2, w2 in terms if e2 != e)) / w, 4)
            if v >= 0:
                new[e].append(v)
        if not new:
            break
        solved.update({e: _mode(vs) for e, vs in new.items()})
    return solved


def fitted_model(db, rules=None, use=None):
    """M5 penalty model (SPEC 5.10): the provisional model plus altnum twins at 0, elements solved from annotations
    with several substitutions, and unobserved substitutions priced at FALLBACK_SCALE x the relation median (an
    analogy the program never chose is evidently expensive). `use` selects the annotations to learn from; default:
    B4-train (annotation_split), so B4-test stays a held-out measure."""
    if use is None:
        train = set(annotation_split(db)["train"])
        use = lambda kind, prm: f"{kind}:{' '.join(prm.types)}" in train     # noqa: E731
    model = provisional_model(db, rules, use)
    for m in MATRICES:
        mat = model.matrix(m)
        for a, b in ALTNUM_TWINS:
            mat.explicit[(a, b)] = mat.explicit[(b, a)] = 0.0
    for (m, a, b), v in solved_elements(db, model, use).items():
        model.matrix(m).explicit[(a, b)] = v
    for m in MATRICES:
        mat = model.matrix(m)
        mat.fallback = {rel: v * FALLBACK_SCALE for rel, v in mat.fallback.items()}
        mat._cache.clear()
    return model

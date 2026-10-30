"""Resolution of resonance (P-I; SPEC 5.3), implemented to reproduce the published trail order.

The algorithm works on the subgraph of AROMATIC input bonds. P-I traverses it depth-first: at a node it takes an
untraversed edge, choosing a bond order where ambiguous; on arriving at a nitrogen it chooses a formal charge (0 or
+1), at a terminal oxygen 0 or -1; a node without untraversed edges has its valence checked; when a branch is complete
it jumps back to the last node that still has untraversed edges ("constructive backtracking"); a valence violation
or a complete structure triggers backtracking to the last open choice ("destructive backtracking"). Each complete
structure is scored and kept only if its score is strictly lower than the best so far; a score of 0 ends the search.

Choice order, read from P-I Fig. 4 (structures H, I, J in that order): the lowest-index untraversed neighbour first,
double bond before single, neutral before charged; the traversal starts at the lowest-index node (SPEC 5.3: to be
calibrated on B5). Because the trail is a stack of choices, the search is a lexicographic enumeration over the choice
sequence, so it is implemented here as a recursion over the fixed sequence of traversal events, pruning any partial
assignment whose valences already fail. Pruning never removes a complete structure, so structures are found in the
same order as by the trail algorithm and the first structure with the lowest score wins.

Score (P-I eq 2): c1 |Z_mol| + c2 sum_{Z<0}|Z| + c3 sum_{Z>0} Z + c4 (N_potential_aromatic - N_aromatic), with
(c1..c4) = (8, 4, 3, 2) since program 0.9.7 and (4, 2, 2, 1) before. The charges scored are those of the subgraph
nodes, i.e. the ones this algorithm chooses (interpretation, FINDINGS 6): counting fixed charges elsewhere in the
molecule (e.g. a phosphate dianion) would make it cheaper to break aromaticity with an N+ that reduces |Z_mol|, which
mistypes nucleotides such as AMP in B2-train.

Not in P-I (fallbacks, flagged with warnings): if no structure is valid, carbon nodes may also be -1 (e.g. an
isoindenide given with aromatic bonds); if still none, or if the search exceeds its step budget, the best structure
found so far is used, or all aromatic bonds become single.
"""
from dataclasses import dataclass, field

from ..molecule import AROMATIC
from .aromaticity import classify_rings

COEFFICIENTS = (8, 4, 3, 2)
COEFFICIENTS_PRE_097 = (4, 2, 2, 1)
MAX_STEPS = 500_000


def allowed_valences(element, charge):
    """Total bond order allowed for an atom of this element and formal charge (None = not checked)."""
    if element == "C":
        return {4} if charge == 0 else {3} if abs(charge) == 1 else set()
    if element == "N":
        return {3 + charge} if -1 <= charge <= 1 else set()
    if element == "O":
        return {2 + charge} if -1 <= charge <= 1 else set()
    if element in ("S", "Se"):
        return {0: {2, 4, 6}, -1: {1}, 1: {3}}.get(charge, set())
    if element == "P":
        return {0: {3, 5}, 1: {4}, -1: {2}}.get(charge, set())
    if element == "B":
        return {0: {3}, -1: {4}}.get(charge, set())
    if element in ("H", "F", "Cl", "Br", "I"):
        return {1} if charge == 0 else set()
    return None


def derived_charge(element, valence):
    """Formal charge implied by an atom's total bond order when the input does not give charges (mol2)."""
    if element == "N":
        return valence - 3
    if element == "O":
        return valence - 2
    if element in ("S", "Se"):
        return {1: -1, 3: 1}.get(valence, 0)
    if element == "C":
        return -1 if valence == 3 else 0
    if element == "B":
        return -1 if valence == 4 else 0
    if element == "P":
        return 1 if valence == 4 else 0
    return 0


@dataclass
class ResonanceResult:
    orders: dict                                  # frozenset((i, j)) -> 1/2/3 for every bond
    charges: list                                 # formal charge of every atom in the chosen structure
    score: float = None                           # None if there was no aromatic subgraph
    structures: list = field(default_factory=list)   # scores of the complete structures, in the order found
    warnings: list = field(default_factory=list)


class _Budget(Exception):
    pass


def _events(nodes, edges, neighbors):
    """Traversal events ('node', i) / ('edge', i, j) for all components, each from its lowest-index node."""
    sub = {i: [j for j in neighbors[i] if frozenset((i, j)) in edges] for i in nodes}
    used, visited, events = set(), set(), []
    for start in sorted(nodes):
        if start in visited:
            continue
        visited.add(start); events.append(("node", start))
        trail, cur = [start], start
        while True:
            nxt = next((j for j in sub[cur] if frozenset((cur, j)) not in used), None)
            if nxt is None:                       # constructive backtracking: last trail node with open edges
                cur = next((t for t in reversed(trail) if any(frozenset((t, j)) not in used for j in sub[t])), None)
                if cur is None:
                    break
                continue
            used.add(frozenset((cur, nxt))); events.append(("edge", cur, nxt))
            if nxt not in visited:
                visited.add(nxt); events.append(("node", nxt))
            trail.append(nxt); cur = nxt
    return events


def resolve(n_atoms, elements, neighbors, bonds, charges, rings, coefficients=COEFFICIENTS, max_steps=MAX_STEPS,
            keep_all=False):
    """Choose bond orders for AROMATIC bonds and formal charges for the nodes of the aromatic subgraph.

    bonds: (i, j, order) as in InputMolecule; charges: input formal charges or None (derived from valences).
    keep_all: also return each complete structure (orders, charges) in ResonanceResult.structures.
    """
    fixed = {frozenset((i, j)): o for i, j, o in bonds if o != AROMATIC}
    edges = {frozenset((i, j)) for i, j, o in bonds if o == AROMATIC}
    nodes = sorted({a for e in edges for a in e})
    base_valence = [sum(o for e, o in fixed.items() if i in e) for i in range(n_atoms)]
    if charges is None:
        node_set = set(nodes)
        charges = [0 if i in node_set else derived_charge(elements[i], base_valence[i]) for i in range(n_atoms)]
    result = ResonanceResult(dict(fixed), list(charges))
    if not edges:
        return result
    for e in edges:
        result.orders[e] = 1

    def options(i, allow_carbanion):
        if elements[i] == "N":
            return (0, 1)
        if elements[i] == "O" and len(neighbors[i]) == 1:
            return (0, -1)
        if elements[i] == "C" and allow_carbanion and charges[i] == 0:
            return (0, -1)
        return (charges[i],)

    events = _events(nodes, edges, neighbors)
    c_mol, c_neg, c_pos, c_arom = coefficients

    def search(allow_carbanion):
        opts = {i: options(i, allow_carbanion) for i in nodes}
        val = list(base_valence)
        remaining = {i: sum(1 for j in neighbors[i] if frozenset((i, j)) in edges) for i in nodes}
        q = dict.fromkeys(nodes)
        orders = {}
        best = {"score": None, "orders": None, "charges": None}
        found = []
        steps = [0]

        def ok(i):
            """Some allowed valence is still reachable: each unassigned edge will add 1 or 2."""
            c = q[i]
            lo, hi = val[i] + remaining[i], val[i] + 2 * remaining[i]
            for x in ([c] if c is not None else opts[i]):
                allowed = allowed_valences(elements[i], x)
                if allowed is None or any(lo <= v <= hi for v in allowed):
                    return True
            return False

        def complete():
            all_q = list(charges)
            for i in nodes:
                all_q[i] = q[i]
            all_orders = {**fixed, **orders}
            classes, n_pot, _ = classify_rings(rings, neighbors, elements,
                                               lambda a, b: all_orders.get(frozenset((a, b)), 0), all_q)
            n_arom = sum(c == "arom" for c in classes.values())
            sub_q = [q[i] for i in nodes]                 # charges outside the subgraph are not scored
            score = (c_mol * abs(sum(sub_q)) + c_neg * sum(-x for x in sub_q if x < 0)
                     + c_pos * sum(x for x in sub_q if x > 0) + c_arom * (n_pot - n_arom))
            found.append((score, dict(orders), dict(q)) if keep_all else score)
            if best["score"] is None or score < best["score"]:
                best.update(score=score, orders=dict(orders), charges=all_q)
            return score == 0

        def step(k):
            steps[0] += 1
            if steps[0] > max_steps:
                raise _Budget
            if k == len(events):
                return complete()
            ev = events[k]
            if ev[0] == "node":
                i = ev[1]
                for c in opts[i]:
                    q[i] = c
                    if ok(i) and step(k + 1):
                        return True
                q[i] = None
                return False
            _, a, b = ev
            e = frozenset((a, b))
            remaining[a] -= 1; remaining[b] -= 1
            try:
                for o in (2, 1):
                    orders[e] = o
                    val[a] += o; val[b] += o
                    stop = ok(a) and ok(b) and step(k + 1)
                    val[a] -= o; val[b] -= o
                    if stop:
                        return True
            finally:
                remaining[a] += 1; remaining[b] += 1
                orders.pop(e, None)
            return False

        truncated = False
        try:
            step(0)
        except _Budget:
            truncated = True
        return best, found, truncated

    best, found, truncated = search(False)
    if truncated:
        result.warnings.append(f"resonance search stopped after {max_steps} steps; best structure so far used")
    if best["score"] is None and not truncated:
        best, found, truncated = search(True)
        if best["score"] is not None:
            result.warnings.append("no valid resonance structure without carbanions; carbon -1 allowed")
    result.structures = found
    if best["score"] is None:
        result.warnings.append("no valid resonance structure: aromatic bonds set to single")
        return result
    result.orders.update(best["orders"])
    result.charges = best["charges"]
    result.score = best["score"]
    return result

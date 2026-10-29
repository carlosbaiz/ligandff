"""Ring perception (P-I; SPEC 5.2).

Rings are cycles of fewer than 8 atoms; larger rings are treated as acyclic. Each non-bridgehead ring atom needs its
smallest ring, each bridgehead atom its three smallest. P-I computes these with a homeomorphically reduced graph plus
a non-SSSR breadth-first message passing; SPEC 5.2 accepts any algorithm giving the same ring sets. Here: the same
first stage (iteratively strip atoms with fewer than two remaining neighbours), then exhaustive enumeration of simple
cycles of 3-7 atoms on what is left, which is complete for that size limit.

Bridgehead (open point, SPEC 5.2): an atom with at least 3 ring bonds, as in the pilot.
"""
MAX_RING = 7


def canonical_ring(atoms):
    """Rotate so the lowest index is first and orient so its lower neighbour comes second."""
    r = list(atoms)
    k = r.index(min(r))
    r = r[k:] + r[:k]
    if r[1] > r[-1]:
        r = [r[0]] + r[1:][::-1]
    return tuple(r)


def ring_bonds(ring):
    return {frozenset((ring[i], ring[(i + 1) % len(ring)])) for i in range(len(ring))}


def _core(neighbors):
    """Atoms that can be in a ring: repeatedly remove atoms with < 2 remaining neighbours (P-I stage 1)."""
    alive = [True] * len(neighbors)
    degree = [len(x) for x in neighbors]
    stack = [i for i, d in enumerate(degree) if d < 2]
    while stack:
        i = stack.pop()
        if not alive[i]:
            continue
        alive[i] = False
        for j in neighbors[i]:
            if alive[j]:
                degree[j] -= 1
                if degree[j] < 2:
                    stack.append(j)
    return alive


def find_rings(neighbors, max_size=MAX_RING):
    """All simple cycles with 3..max_size atoms, as canonical tuples, sorted by (size, atoms)."""
    alive = _core(neighbors)
    found = set()

    def extend(start, path, on_path):
        cur = path[-1]
        for nb in neighbors[cur]:
            if nb == start and len(path) >= 3:
                found.add(canonical_ring(path))
            elif nb > start and alive[nb] and nb not in on_path and len(path) < max_size:
                path.append(nb); on_path.add(nb)
                extend(start, path, on_path)
                path.pop(); on_path.discard(nb)

    for s in range(len(neighbors)):
        if alive[s]:
            extend(s, [s], {s})
    return sorted(found, key=lambda r: (len(r), r))


def atom_rings(n_atoms, neighbors, rings):
    """(per-atom ring lists, bridgehead flags): smallest ring, or three smallest for bridgeheads; ties by atoms."""
    rbonds = set().union(*map(ring_bonds, rings)) if rings else set()
    bridgehead = [sum(frozenset((i, j)) in rbonds for j in neighbors[i]) >= 3 for i in range(n_atoms)]
    member = [[] for _ in range(n_atoms)]
    for r in rings:                                         # rings are sorted by (size, atoms)
        for a in r:
            member[a].append(r)
    return [m[:3] if bridgehead[i] else m[:1] for i, m in enumerate(member)], bridgehead

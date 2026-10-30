"""Ring classes and aromaticity (P-I; SPEC 5.4).

Classes: 'arom'; 'sp2' (every ring atom in a double bond, except at most one N/O/P/S with only single bonds);
'sp3' (no ring atom in any double bond); 'mixed' (everything else). Aromatic = 5-7-membered ring whose atoms all
have <= 3 substituents ("potential aromatic") with exactly 6 in-ring pi electrons:
  in-ring double/triple bond 2 (1 per atom); atom with an exocyclic double bond that is in another aromatic ring 1;
  heteroatom with only single bonds 2, or 1 or 2 (whichever makes the ring aromatic) if in another aromatic ring.
Iterated over all potential aromatic rings until no ring changes status.
Extension kept from the pilot (calibration point): a carbanion with only single bonds counts 2 (cyclopentadienide).
"""
MAX_ITERATIONS = 50


def _is_potential(ring, neighbors):
    return 5 <= len(ring) <= 7 and all(len(neighbors[a]) <= 3 for a in ring)


def classify_rings(rings, neighbors, elements, order, charges):
    """Return ({ring: class}, number of potential aromatic rings, converged).

    order(i, j) -> resolved bond order (1, 2, 3); charges: resonance formal charges.
    """
    potential = [r for r in rings if _is_potential(r, neighbors)]
    in_ring = {r: {frozenset((r[k], r[(k + 1) % len(r)])) for k in range(len(r))} for r in potential}
    aromatic, converged = set(), False
    for _ in range(MAX_ITERATIONS):
        new = set()
        for r in potential:
            fixed = flexible = 0
            for a in r:
                multiple = [b for b in neighbors[a] if order(a, b) >= 2]
                in_other = any(a in r2 for r2 in aromatic if r2 != r)
                if any(frozenset((a, b)) in in_ring[r] for b in multiple):
                    fixed += 1
                elif multiple:
                    fixed += 1 if in_other else 0
                elif elements[a] not in ("C", "B"):
                    if in_other:
                        fixed += 1; flexible += 1
                    else:
                        fixed += 2
                elif elements[a] == "C" and charges[a] < 0:
                    fixed += 2
            if fixed <= 6 <= fixed + flexible:
                new.add(r)
        if new == aromatic:
            converged = True
            break
        aromatic = new
    classes = {}
    for r in rings:
        if r in aromatic:
            classes[r] = "arom"
            continue
        double = [any(order(a, b) >= 2 for b in neighbors[a]) for a in r]
        if not any(double):
            classes[r] = "sp3"
        elif sum(not d for d in double) <= 1 and all(d or elements[a] in ("N", "O", "P", "S") for a, d in zip(r, double)):
            classes[r] = "sp2"
        else:
            classes[r] = "mixed"
    return classes, len(potential), converged

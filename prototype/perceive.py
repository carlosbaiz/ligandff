"""CGenFF-style perception for the pilot typer (Paper I): rings <= 7, lowest-penalty resonance structure
(8|Zmol| + 4 sum|Z-| + 3 sum|Z+| + 2(N_potential_aromatic - N_aromatic)), ring classes
(aromatic by Paper I's electron counting, all-sp2, all-sp3, mixed)."""
from rdkit import Chem

BO = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}


class Mol:
    pass


def _cycles(nbrs, heavy):
    found = {}
    def dfs(start, cur, path):
        for nb in nbrs[cur]:
            if nb == start and len(path) >= 3:
                r = path[:]
                i = r.index(min(r)); r = r[i:] + r[:i]
                if r[1] > r[-1]: r = [r[0]] + r[1:][::-1]
                found[tuple(r)] = True
            elif nb > start and nb not in path and len(path) < 7 and heavy[nb]:
                path.append(nb); dfs(start, nb, path); path.pop()
    for s in range(len(nbrs)):
        if heavy[s]: dfs(s, s, [s])
    return list(found)


def _ring_bonds(r):
    return {frozenset((r[i], r[(i + 1) % len(r)])) for i in range(len(r))}


def _classify(m, bo, charge):
    """Return dict ring -> class ('arom', 'sp2', 'sp3', 'mixed') following Paper I."""
    rings = m.rings
    pot = [r for r in rings if 5 <= len(r) <= 7 and all(m.deg[a] <= 3 for a in r)]
    arom = set()
    for _ in range(10):
        new = set()
        for r in pot:
            rb = _ring_bonds(r)
            fixed, flex = 0, 0
            for a in r:
                in_ring_db = any(bo[a][b] >= 2 and frozenset((a, b)) in rb for b in m.nbrs[a])
                out_db = any(bo[a][b] >= 2 and frozenset((a, b)) not in rb for b in m.nbrs[a])
                other_arom = any(a in r2 for r2 in arom if r2 != r)
                if in_ring_db: fixed += 1
                elif out_db: fixed += 1 if other_arom else 0
                elif m.el[a] != 'C' and m.el[a] != 'B':
                    if other_arom: fixed += 1; flex += 1
                    else: fixed += 2
                elif m.el[a] == 'C' and charge[a] < 0: fixed += 2
            if fixed <= 6 <= fixed + flex: new.add(r)
        if new == arom: break
        arom = new
    cls = {}
    for r in rings:
        if r in arom: cls[r] = 'arom'; continue
        dbl = [any(bo[a][b] >= 2 for b in m.nbrs[a]) for a in r]
        if not any(dbl): cls[r] = 'sp3'
        elif sum(not d for d in dbl) <= 1 and all(d or m.el[a] in 'NOPS' for a, d in zip(r, dbl)): cls[r] = 'sp2'
        else: cls[r] = 'mixed'
    return cls, len(pot)


def perceive(rdmol, max_res=64):
    m = Mol()
    m.n = rdmol.GetNumAtoms()
    m.el = [a.GetSymbol() for a in rdmol.GetAtoms()]
    m.nbrs = [[x.GetIdx() for x in a.GetNeighbors()] for a in rdmol.GetAtoms()]
    m.deg = [len(x) for x in m.nbrs]
    m.nH = [sum(m.el[j] == 'H' for j in m.nbrs[i]) for i in range(m.n)]
    heavy = [e != 'H' for e in m.el]
    m.rings = _cycles(m.nbrs, heavy)
    # resonance structures (Kekule forms); pick lowest CGenFF penalty
    cands = []
    try:
        sup = Chem.ResonanceMolSupplier(rdmol, Chem.KEKULE_ALL)
        for k, rm in enumerate(sup):
            if rm is None or k >= max_res: continue
            cands.append(rm)
    except Exception:
        pass
    if not cands:
        k = Chem.Mol(rdmol); Chem.Kekulize(k, clearAromaticFlags=True); cands = [k]
    best = None
    for rm in cands:
        bo = [dict() for _ in range(m.n)]
        for b in rm.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            o = BO.get(b.GetBondType(), 1)
            if b.GetBondType() == Chem.BondType.AROMATIC: o = 1
            bo[i][j] = bo[j][i] = o
        charge = [a.GetFormalCharge() for a in rm.GetAtoms()]
        cls, npot = _classify(m, bo, charge)
        narom = sum(v == 'arom' for v in cls.values())
        ringpos = sum(1 for a in range(m.n) if charge[a] > 0 and any(a in r for r in m.rings))
        pen = -0.001 * ringpos + 8 * abs(sum(charge)) + 4 * sum(-c for c in charge if c < 0) + 3 * sum(c for c in charge if c > 0) + 2 * (npot - narom)
        if best is None or pen < best[0]: best = (pen, bo, charge, cls)
    _, m.bo, m.charge, m.rcls = best
    m.val = [sum(m.bo[i].values()) for i in range(m.n)]
    # per-atom ring lists: smallest ring for ordinary ring atoms, three smallest for bridgeheads (>=3 ring bonds)
    m.atom_rings = [[] for _ in range(m.n)]
    ringbonds = set().union(*[_ring_bonds(r) for r in m.rings]) if m.rings else set()
    for i in range(m.n):
        mine = sorted([r for r in m.rings if i in r], key=len)
        nrb = sum(frozenset((i, j)) in ringbonds for j in m.nbrs[i])
        m.atom_rings[i] = mine[:3] if nrb >= 3 else mine[:1]
    m.inring = lambda i, j: frozenset((i, j)) in ringbonds
    # atoms that carry a positive / negative formal charge in at least one resonance structure
    m.posany = [c > 0 for c in m.charge]; m.negany = [c < 0 for c in m.charge]
    try:
        for k, rm in enumerate(Chem.ResonanceMolSupplier(rdmol)):
            if rm is None or k >= max_res: continue
            for a in rm.GetAtoms():
                if a.GetFormalCharge() > 0: m.posany[a.GetIdx()] = True
                if a.GetFormalCharge() < 0: m.negany[a.GetIdx()] = True
    except Exception:
        pass
    # fused/bridged ring pairs (share >= 2 atoms) per atom
    m.fused = [any(len(set(r1) & set(r2)) >= 2 and sum(j in set(r1) | set(r2) for j in m.nbrs[i]) >= 3
                   for a, r1 in enumerate(m.atom_rings[i]) for r2 in m.atom_rings[i][a + 1:]) for i in range(m.n)]
    return m

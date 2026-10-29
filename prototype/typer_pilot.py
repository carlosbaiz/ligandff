"""Pilot CGenFF-compatible atom typer (version 4 of the pilot; see docs/FINDINGS.md).
Python rules + RDKit, written from the MASS-line type descriptions, Paper I, and three rounds of
training-split failure analysis. Reference only: the production typer should use declarative rules
(docs/SPEC.md, M3). Entry point: type_molecule(rdkit_mol_with_explicit_H) -> list of CGenFF types."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perceive import perceive

HET = set('N O S P'.split())


def type_molecule(rdmol):
    m = perceive(rdmol)
    T = [None] * m.n
    ctx = Ctx(m)
    for i in range(m.n):
        if m.el[i] != 'H': T[i] = ctx.heavy(i)
    for i in range(m.n):
        if m.el[i] == 'H': T[i] = ctx.hydrogen(i, T)
    return T


class Ctx:
    def __init__(s, m):
        s.m = m

    # ---------- helpers ----------
    def rings(s, i, cls=None, size=None):
        return [r for r in s.m.atom_rings[i] if (cls is None or s.m.rcls[r] == cls) and (size is None or len(r) == size)]

    def arom(s, i, size=None): return bool(s.rings(i, 'arom', size))
    def planar5(s, i): return [r for r in s.m.atom_rings[i] if len(r) == 5 and s.m.rcls[r] in ('arom', 'sp2')]
    def heavy_nb(s, i): return [j for j in s.m.nbrs[i] if s.m.el[j] != 'H']
    def dbl(s, i): return [j for j in s.m.nbrs[i] if s.m.bo[i][j] == 2]
    def tpl(s, i): return [j for j in s.m.nbrs[i] if s.m.bo[i][j] == 3]
    def is_carbonyl_c(s, i): return s.m.el[i] == 'C' and any(s.m.el[j] in ('O', 'S') and s.m.bo[i][j] == 2 for j in s.m.nbrs[i])
    def cat_amidine_c(s, k):
        m = s.m
        nN = [x for x in m.nbrs[k] if m.el[x] == 'N' and not s.nitro_like(x)]
        return m.el[k] == 'C' and len(nN) >= 2 and any(m.bo[k][x] == 2 for x in nN) and any(m.charge[x] > 0 for x in nN) \
            and not s.arom(k)
    def posN(s, j):
        m = s.m
        if m.el[j] != 'N': return False
        if m.charge[j] > 0 and not any(m.charge[b] < 0 for b in m.nbrs[j]): return True
        if any(s.ring_cationic(r) for r in m.atom_rings[j] if s.m.rcls[r] in ('arom', 'sp2')): return True
        return any(s.cat_amidine_c(k) for k in m.nbrs[j])
    def nF(s, i): return sum(s.m.el[j] == 'F' for j in s.m.nbrs[i])
    def nH(s, i): return s.m.nH[i]
    def same_ring(s, i, j, size=None):
        return any(i in r and j in r and (size is None or len(r) == size) for r in s.m.rings)
    def ring_cationic(s, r): return any(s.m.charge[a] > 0 and not any(s.m.charge[b] < 0 for b in s.m.nbrs[a]) for a in r)
    def link5(s, i): return [k for k in s.heavy_nb(i) if s.m.bo[i][k] == 1 and not s.same_ring(i, k)
                             and s.m.atom_rings[k] and min(len(r) for r in s.m.atom_rings[k]) == 5]
    def nitro_like(s, k): return s.m.el[k] == 'N' and sum(s.m.el[x] == 'O' for x in s.m.nbrs[k]) >= 2

    # ---------- heavy atoms ----------
    def heavy(s, i):
        e = s.m.el[i]
        return {'C': s.carbon, 'N': s.nitrogen, 'O': s.oxygen, 'S': s.sulfur, 'P': s.phosphorus, 'F': s.halogen,
                'Cl': s.halogen, 'Br': s.halogen, 'I': s.halogen, 'B': lambda i: 'BG201', 'Se': lambda i: 'SEGD1'}.get(e, lambda i: 'UNK')(i)

    def carbon(s, i):
        m = s.m
        if s.tpl(i):
            j = s.tpl(i)[0]
            if m.el[j] == 'N': return 'CG1N1'
            return 'CG1T2' if s.nH(i) else 'CG1T1'
        d = s.dbl(i)
        hv = s.heavy_nb(i)
        if d:
            j = d[0]
            if m.el[j] in ('O', 'S') and not s.same_ring(i, j):     # carbonyl / thiocarbonyl
                if m.el[j] == 'O' and sum(m.el[k] == 'O' and m.bo[i][k] == 2 for k in m.nbrs[i]) == 2: return 'CG2O7'
                if s.arom(i, 6): return 'CG2R63'
                if any(m.el[k] in ('N', 'O', 'S') and k != j and (s.same_ring(i, k, 5) or s.same_ring(i, k, 4)) for k in hv):
                    return 'CG2R53'
                if s.rings(i, None, 4) and m.el[j] == 'O': return 'CG2O4'
                singles = [k for k in hv if k != j]
                nN = sum(m.el[k] == 'N' for k in singles)
                nO = sum(m.el[k] in ('O', 'S') for k in singles)
                if any(m.el[k] == 'O' and m.negany[k] for k in singles): return 'CG2O3'
                if nN + nO >= 2: return 'CG2O6'
                if nN == 1: return 'CG2O1'
                if nO == 1: return 'CG2O2'
                if s.nH(i): return 'CG2O4'
                return 'CG2O5'
            # six-membered aromatic carbons
            if s.arom(i, 6):
                if s.nF(i): return 'CG2R66'
                if any(m.bo[i][k] == 1 and s.arom(k, 6) and m.el[k] == 'C' and
                       not any(i in r and k in r for r in s.rings(i, 'arom', 6)) for k in hv): return 'CG2R67'
                if s.rings(i, None, 5) or s.rings(i, None, 4): return 'CG2RC0'
                nN = [k for k in hv if m.el[k] == 'N' and not s.nitro_like(k)]
                if len(nN) >= 2: return 'CG2R64'
                r6 = s.rings(i, 'arom', 6)
                if any(s.ring_cationic(r) or any(s.is_carbonyl_c(a) for a in r) for r in r6): return 'CG2R62'
                return 'CG2R61'
            # five-membered planar rings (aromatic, all-sp2) or C=X double bond inside a 5-ring
            p5 = s.planar5(i)
            in5 = any(s.same_ring(i, k, 5) for k in d)
            if s.rings(i, None, 4) and any(s.same_ring(i, k, 4) for k in d): return 'CG2RC0'
            if p5 or in5:
                if s.arom(i, 6): return 'CG2RC0'
                exo = not any(s.same_ring(i, k, 5) for k in d) and not s.arom(i, 5)
                if exo and m.el[j] == 'C' and not s.arom(j):
                    return 'CG251O' if any(m.el[k] in HET for k in hv) else 'CG25C1'
                if s.link5(i): return 'CG2R57'
                dN = [k for k in d if m.el[k] == 'N']
                if dN:
                    other = [k for k in hv if m.el[k] in HET and k not in dN]
                    return 'CG2R53' if other else 'CG2R52'
                return 'CG2R51'
            if s.arom(i, 7): return 'CG2R71'
            # amidine / guanidine
            nNnb = [k for k in hv if m.el[k] == 'N' and not s.nitro_like(k)]
            if m.el[j] == 'N' and len(nNnb) == 3: return 'CG2N1'
            if m.el[j] == 'N' and len(nNnb) == 2 and any(m.posany[k] for k in nNnb): return 'CG2N2'
            if s.rings(i, None, 5) and not s.same_ring(i, j):
                return 'CG251O' if any(m.el[k] in HET for k in hv) else 'CG25C1'
            def conj_from(a, skip):
                return any(m.bo[a][k] == 1 and m.el[k] != 'H' and any(m.bo[k][l] >= 2 for l in m.nbrs[k] if l != a) for k in m.nbrs[a] if k != skip)
            conj = conj_from(i, j) or conj_from(j, i)
            het = any(m.el[k] in HET for k in hv if k != j)
            if conj:
                if s.nH(i) == 2: return 'CG2DC3'
                return 'CG2D1O' if het else 'CG2DC1'
            if het: return 'CG2D1O'
            return 'CG2D2' if s.nH(i) == 2 else 'CG2D1'
        # aromatic 5-ring carbon whose Kekule double bond lies outside the ring (fusion atoms)
        if s.arom(i, 5) or s.arom(i, 6):
            return 'CG2RC0'
        # sp3 carbon
        if m.fused[i] and any(len(r) <= 5 for r in m.atom_rings[i]): return 'CG3RC1'
        if s.rings(i, None, 3): return 'CG3C31'
        if s.rings(i, None, 4): return 'CG3C41'
        f = s.nF(i)
        posN = any(s.posN(k) for k in m.nbrs[i])
        h = s.nH(i)
        if s.rings(i, None, 5):
            if h == 0: return 'CG3C50'
            if h == 1: return 'CG3C53' if posN else 'CG3C51'
            return 'CG3C54' if posN else 'CG3C52'
        if f == 1: return 'CG322'
        if f == 2: return 'CG312'
        if f >= 3: return 'CG302'
        if h == 0: return 'CG301'
        if h == 1: return 'CG314' if posN else 'CG311'
        if h == 2:
            if any(m.el[k] == 'S' and m.charge[k] < 0 for k in m.nbrs[i]): return 'CG323'
            return 'CG324' if posN else 'CG321'
        return 'CG334' if posN else 'CG331'

    def nitrogen(s, i):
        m = s.m
        hv = s.heavy_nb(i)
        if s.tpl(i): return 'NG1T1'
        if s.nitro_like(i): return 'NG2O1'
        if m.charge[i] > 0 and m.val[i] == 4 and not s.dbl(i) and not s.arom(i):
            return {0: 'NG3P0', 1: 'NG3P1', 2: 'NG3P2', 3: 'NG3P3'}[min(s.nH(i), 3)]
        d = s.dbl(i)
        p5 = s.planar5(i)
        if p5 and not s.arom(i, 6):
            if m.posany[i] and any(s.ring_cationic(r) for r in p5): return 'NG2R52'
            if d: return 'NG2R50'
            if any(s.arom(k, 5) and not m.inring(i, k) for k in hv): return 'NG2R57'
            return 'NG2R51'
        if s.arom(i, 6):
            if s.planar5(i): return 'NG2RC0'
            if d and m.charge[i] <= 0:
                near = set()
                for k in m.nbrs[i]:
                    near.add(k); near.update(m.nbrs[k])
                near.discard(i)
                return 'NG2R62' if any(m.el[a] in HET and s.arom(a) for a in near) else 'NG2R60'
            if any(s.arom(k, 6) and not m.inring(i, k) for k in hv): return 'NG2R67'
            return 'NG2R61'
        if d:
            if s.posN(i): return 'NG2R52' if s.rings(i, None, 5) else 'NG2P1'
            if any(s.same_ring(i, k, 5) for k in d): return 'NG2R50'
            return 'NG2D1'
        # single-bonded N
        if any(s.cat_amidine_c(k) for k in hv): return 'NG2R52' if s.rings(i, None, 5) else 'NG2P1'
        amideC = [k for k in hv if s.is_carbonyl_c(k)]
        if amideC:
            if any(s.same_ring(i, k, 4) for k in amideC): return 'NG2R43'
            if any(s.same_ring(i, k, 5) for k in amideC): return 'NG2R53'
            return {0: 'NG2S0', 1: 'NG2S1'}.get(s.nH(i), 'NG2S2')
        if s.nH(i) == 2 and len(hv) == 1 and s.arom(hv[0], 6): return 'NG2S3'
        if s.rings(i, None, 5): return 'NG3C51'
        if any(m.el[k] == 'N' for k in hv): return 'NG3N1'
        return {0: 'NG301', 1: 'NG311', 2: 'NG321'}.get(s.nH(i), 'NG331')

    def oxygen(s, i):
        m = s.m
        hv = s.heavy_nb(i)
        if not hv: return 'OGTIP3'
        j = hv[0]
        if m.bo[i][j] == 2 or (len(m.nbrs[i]) == 1 and m.charge[i] < 0):
            if m.el[j] == 'C':
                ct = s.carbon(j)
                if ct == 'CG2O4' and not m.nH[j]: return 'OG2D3'
                return {'CG2O3': 'OG2D2', 'CG2O5': 'OG2D3', 'CG2R63': 'OG2D4', 'CG2O7': 'OG2D5', 'CG2R53': 'OG2D1'}.get(ct, 'OG2D1' if m.bo[i][j] == 2 else 'OG312')
            if m.el[j] == 'N': return 'OG2N1'
            if m.el[j] in ('P', 'S'): return 'OG2P1'
            return 'OG2D1'
        if len(hv) == 1:
            return 'OG311'
        if s.arom(i, 5): return 'OG2R50'
        if s.rings(i, None, 3): return 'OG3C31'
        if s.rings(i, None, 5) or s.rings(i, None, 4): return 'OG3C51'
        if any(m.el[k] in ('P', 'S') for k in hv):
            return 'OG304' if all(m.el[k] in ('P', 'S') for k in hv) else 'OG303'
        if any(s.is_carbonyl_c(k) for k in hv): return 'OG3R60' if s.rings(i, None, 6) else 'OG302'
        if s.rings(i, None, 6):
            return 'OG3R60' if any(s.dbl(k) or s.arom(k) for k in hv) else 'OG3C61'
        return 'OG301'

    def sulfur(s, i):
        m = s.m
        hv = s.heavy_nb(i)
        if s.arom(i, 5): return 'SG2R50'
        nO = sum(m.el[k] == 'O' for k in hv)
        if nO >= 3: return 'SG3O1'
        if nO == 2: return 'SG3O2'
        if nO == 1 and len(hv) >= 3: return 'SG3O3'
        if any(m.bo[i][k] == 2 and m.el[k] == 'C' for k in hv): return 'SG2D1'
        if m.charge[i] < 0: return 'SG302'
        if any(m.el[k] == 'S' for k in hv): return 'SG301'
        return 'SG311'

    def phosphorus(s, i):
        neg = sum(s.m.charge[k] < 0 for k in s.m.nbrs[i])
        return {0: 'PG0', 1: 'PG1'}.get(neg, 'PG2')

    def halogen(s, i):
        m = s.m
        e = m.el[i]
        j = s.heavy_nb(i)[0]
        if s.arom(j): return {'F': 'FGR1', 'Cl': 'CLGR1', 'Br': 'BRGR1', 'I': 'IGR1'}[e]
        n = sum(m.el[k] == e for k in m.nbrs[j])
        if e == 'F': return {1: 'FGA1', 2: 'FGA2'}.get(n, 'FGA3')
        if e == 'Cl': return 'CLGA3' if n >= 3 else 'CLGA1'
        if e == 'Br': return {1: 'BRGA1', 2: 'BRGA2'}.get(n, 'BRGA3')
        return 'IGR1'

    # ---------- hydrogens ----------
    def hydrogen(s, i, T):
        m = s.m
        j = m.nbrs[i][0]
        e, tj = m.el[j], T[j]
        if e == 'C':
            if tj.startswith('CG2R6') or (tj == 'CG2RC0' and s.arom(j, 6)):
                rs = s.rings(j, 'arom', 6)
                if any(s.ring_cationic(r) for r in rs): return 'HGR63'
                if tj in ('CG2R62', 'CG2R64'): return 'HGR62'
                if any(m.el[k] in HET or T[k] == 'CG2R66' for k in s.heavy_nb(j)): return 'HGR62'
                return 'HGR61'
            if tj == 'CG2RC0' and s.rings(j, None, 4) and not s.arom(j): return 'HGA4'
            if tj.startswith('CG2R5') or tj == 'CG2RC0':
                if any(m.el[k] in HET for k in s.heavy_nb(j)) or any(s.ring_cationic(r) for r in s.planar5(j)): return 'HGR52'
                return 'HGR51'
            if tj.startswith('CG2R7'): return 'HGR71'
            if tj in ('CG2O4', 'CG2O1'): return 'HGR52'
            if tj == 'CG1T2': return 'HGPAM1'
            if tj.startswith('CG2'): return 'HGA5' if m.nH[j] == 2 else 'HGA4'
            if any(T[k] == 'NG3P0' for k in m.nbrs[j]): return 'HGP5'
            f = s.nF(j)
            if f == 1: return 'HGA6'
            if f == 2: return 'HGA7'
            return {1: 'HGA1', 2: 'HGA2'}.get(m.nH[j], 'HGA3')
        if e == 'N':
            if s.posN(j) or (m.charge[j] > 0): return 'HGP2'
            if tj == 'NG2S3': return 'HGP4'
            if any(m.el[k] in ('S', 'O', 'N') for k in s.heavy_nb(j)): return 'HGP1'
            if tj == 'NG311': return 'HGPAM1'
            if tj == 'NG321': return 'HGPAM2'
            if tj == 'NG331': return 'HGPAM3'
            return 'HGP1'
        if e == 'S': return 'HGP3'
        return 'HGP1'

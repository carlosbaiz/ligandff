"""How much of the CGenFF substitution penalty matrix can be recovered from the
'from <source>, penalty= N' annotations in the public parameter file?

Paper II: TPS = sum of per-position substitution penalties; inner positions (bond atoms, angle centre,
dihedral 2/3, improper centre) use the 'bonded' matrix x10; outer positions use the 'nonbonded' matrix
(impropers use the bonded matrix for all four, x10 only for the centre); plus bond-group penalties.
For annotations with exactly one substituted position, penalty/scale is one matrix element (if no bond-group
term applied). Bond-group terms are >= 20, so small penalties are clean."""
import re, sys, collections

prm = sys.argv[1]
N = {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}
sec = None
obs = collections.defaultdict(list)   # (matrix, from_type, to_type) -> list of element values
nsingle = collections.Counter(); nall = 0; multi = 0
seen = set()
for l in open(prm, errors='ignore'):
    s = l.split('!')[0].split()
    if s and s[0].upper() in ('BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS', 'IMPROPER', 'NONBONDED', 'NBFIX', 'END'):
        sec = s[0].upper().replace('IMPROPER', 'IMPROPERS').replace('IMPROPERSS', 'IMPROPERS'); continue
    if sec not in N: continue
    m = re.search(r'from ([A-Z0-9 ]+?),\s*penalty=\s*([\d.]+)', l)
    if not m: continue
    tgt = s[:N[sec]]; src = m.group(1).split(); p = float(m.group(2))
    if len(src) != len(tgt): continue
    key = (sec, tuple(tgt), tuple(src))
    if key in seen: continue
    seen.add(key); nall += 1
    best = None
    for cand in (src, src[::-1]):          # the program also tries the reversed source
        diffs = [i for i, (a, b) in enumerate(zip(tgt, cand)) if a != b]
        if best is None or len(diffs) < len(best[0]): best = (diffs, cand)
    diffs, cand = best
    if len(diffs) != 1: multi += 1; continue
    i = diffs[0]; n = N[sec]
    if sec == 'BONDS': mat, scale = 'bonded', 10
    elif sec == 'IMPROPERS': mat, scale = 'bonded', (10 if i == 0 else 1)
    elif sec == 'ANGLES': mat, scale = ('bonded', 10) if i == 1 else ('nonbonded', 1)
    else: mat, scale = ('bonded', 10) if i in (1, 2) else ('nonbonded', 1)
    nsingle[sec] += 1
    obs[(mat, tgt[i], cand[i])].append(p / scale)

print(f"distinct annotated parameters: {nall}; single-substitution: {sum(nsingle.values())} {dict(nsingle)}; multi-substitution: {multi}")
types = set()
for (mat, a, b) in obs: types |= {a, b}
for mat in ('bonded', 'nonbonded'):
    els = {k: v for k, v in obs.items() if k[0] == mat}
    vals = [min(v) for v in els.values()]   # min: other occurrences may include bond-group terms
    consistent = sum(1 for v in els.values() if max(v) - min(v) < 1e-6)
    print(f"{mat:9s} matrix: {len(els)} distinct off-diagonal elements observed "
          f"({consistent} always the same value; {len(els)-consistent} vary -> bond groups or version changes)")
print(f"atom types involved: {len(types)}  (full matrices would have ~{166*165} off-diagonal elements each)")

# check against the NG3 extract published in Paper II Table 1 (one of the two matrices; which one is not stated)
tab = {('NG3P2','NG3P1'):1,('NG3P2','NG3P3'):3,('NG3P2','NG3P0'):4,('NG3P3','NG3P2'):1,('NG3P3','NG3P1'):2,('NG3P3','NG3P0'):4,
       ('NG3P1','NG3P2'):1,('NG3P1','NG3P0'):3,('NG3P1','NG3P3'):4,('NG3P0','NG3P1'):1,('NG3P0','NG3P2'):2,('NG3P0','NG3P3'):4,
       ('NG321','NG311'):1,('NG321','NG301'):1.5,('NG321','NG3N1'):2.5,('NG321','NG3C51'):3,('NG321','NG331'):4,
       ('NG311','NG301'):.5,('NG311','NG321'):1,('NG311','NG3N1'):1.5,('NG311','NG3C51'):2,('NG311','NG331'):5,
       ('NG301','NG311'):.5,('NG301','NG321'):1.5,('NG301','NG3N1'):2,('NG301','NG3C51'):2.5,('NG301','NG331'):5.5,
       ('NG3C51','NG311'):2,('NG3C51','NG301'):2.5,('NG3C51','NG321'):3,('NG3C51','NG331'):7}
print("\nPaper II Table 1 (NG3 extract) vs annotation-derived elements (from -> to):")
for (a, b), v in tab.items():
    for mat in ('bonded', 'nonbonded'):
        if (mat, a, b) in obs:
            print(f"   {a}->{b} [{mat}] table {v}  observed {sorted(set(obs[(mat,a,b)]))}")

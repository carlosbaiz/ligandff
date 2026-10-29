"""Prototype CGenFF-compatible charge assignment (bond/angle/dihedral charge increments).

Implements Vanommeslaeghe, Raman & MacKerell, JCIM 2012, 52, 3155 (Paper II):
    q_i = q0_i - sum_j beta_ij, extended with two increments per angle and three per dihedral.
    For a parameter instance A-B(-C(-D)) with increments (d1[, d2[, d3]]) in parameter-file orientation,
    each increment is a transfer along one bond: d1 on A-B, d2 on B-C, d3 on C-D, with
    q_first -= d, q_second += d. If the instance matches the parameter in reverse order, the
    order of the atoms is reversed (equivalently: order and sign of the increments reversed).
Per-charge penalty, Paper II eq 4:
    sigma_i = sqrt( sum_n cbrt(|beta_n| + delta) * sigma_n^2 ),  delta = 0.05**6
    where, as verified against CGenFF program 5.0 output, beta_n is the NET charge that one
    parameter instance transfers onto atom i and sigma_n is that parameter's analogy penalty.

Increments come from the public CGenFF 5.0 file par_all36_cgenff_increments.prm (the !<% ... %>! field;
'N/A' = increment not used). This prototype does NOT choose analog parameters itself: missing
parameters must be supplied as {missing types -> (source types, penalty)}, e.g. from a CGenFF .str file.

usage:
  python charges.py fixtures --toppar DIR --fixtures tests/fixtures/si_v5_reference.json
  python charges.py training --toppar DIR      # formula check on all training residues (formal charges inferred)
"""
import argparse, collections, itertools, json, math, re, statistics, sys

KINDS = {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4}
DELTA = 0.05 ** 6
MISSING = object()


def load_increments(prm_path):
    """Return {n_atoms: {types_tuple: [increments] or None}} from a CGenFF increments parameter file."""
    table = {2: {}, 3: {}, 4: {}}
    sec = None
    for line in open(prm_path, errors='ignore'):
        words = line.split('!')[0].split()
        if words and words[0].upper() in ('BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'IMPROPERS', 'NONBONDED', 'NBFIX', 'HBOND', 'END', 'CMAP'):
            sec = words[0].upper(); continue
        if sec not in KINDS or not words: continue
        m = re.search(r'<%(.*?)%>', line)
        if not m: continue
        key = tuple(words[:KINDS[sec]])
        if key in table[KINDS[sec]]: continue                    # multi-term dihedrals repeat the same increments
        body = m.group(1)
        table[KINDS[sec]][key] = None if 'N/A' in body else [float(x) for x in body.split()]
    return table


def _lookup(table, key):
    t = table[len(key)]
    if key in t: return t[key], 1
    r = tuple(reversed(key))
    if r in t: return t[r], -1
    return MISSING, 0


def _instances(n, bonds):
    adj = [set() for _ in range(n)]
    for i, j in bonds: adj[i].add(j); adj[j].add(i)
    out = [tuple(b) for b in bonds]
    for b in range(n):
        for a, c in itertools.combinations(sorted(adj[b]), 2): out.append((a, b, c))
    for b, c in bonds:
        for a in adj[b] - {c}:
            for d in adj[c] - {b, a}: out.append((a, b, c, d))
    return out


def assign(types, bonds, q0, table, analogs=None, max_dihedral_penalty=50.0):
    """types: atom types; bonds: (i, j) index pairs; q0: formal charges; analogs: {missing_types: (source_types, penalty)}.
    Returns (charges, charge_penalties, missing_parameters)."""
    n = len(types)
    analogs = {tuple(k): (tuple(v[0]), float(v[1])) for k, v in (analogs or {}).items()}
    q = [float(x) for x in q0]
    contrib = collections.defaultdict(list)                      # atom -> [(net transfer, penalty)]
    missing = set()
    for inst in _instances(n, bonds):
        key = tuple(types[a] for a in inst)
        vals, orient = _lookup(table, key)
        pen = 0.0
        if vals is MISSING:
            src = analogs.get(key) or analogs.get(tuple(reversed(key)))
            if src is None: missing.add(key); continue
            o1 = 1 if key in analogs else -1                     # instance orientation vs the missing key as written
            vals, o2 = _lookup(table, src[0])
            if vals is MISSING: missing.add(key); continue
            orient, pen = o1 * o2, src[1]
            if len(inst) == 4 and pen > max_dihedral_penalty:    # Paper II: ignore high-penalty dihedral increments
                vals, pen = [0.0, 0.0, 0.0], max_dihedral_penalty
        if vals is None: continue                                 # N/A increments
        atoms = inst if orient == 1 else tuple(reversed(inst))
        net = collections.defaultdict(float)
        for k, d in enumerate(vals):
            net[atoms[k]] -= d; net[atoms[k + 1]] += d
        for a, dq in net.items():
            q[a] += dq
            if pen > 0: contrib[a].append((dq, pen))
    pens = [math.sqrt(sum((abs(dq) + DELTA) ** (1 / 3) * p * p for dq, p in contrib[i])) for i in range(n)]
    return q, pens, sorted(missing)


def cmd_fixtures(args):
    table = load_increments(f"{args.toppar}/par_all36_cgenff_increments.prm")
    ok_all = True
    for mol in json.load(open(args.fixtures)):
        names = [a['name'] for a in mol['atoms']]; idx = {n: i for i, n in enumerate(names)}
        types = [a['type'] for a in mol['atoms']]
        bonds = [(idx[a], idx[b]) for a, b in mol['bonds']]
        q0 = [mol['formal_charges'].get(n, 0) for n in names]
        analogs = {tuple(x['missing']): (x['source'], x['penalty']) for x in mol['analog_parameters']}
        q, pen, missing = assign(types, bonds, q0, table, analogs)
        dq = max(abs(round(a, 3) - b['charge']) for a, b in zip(q, mol['atoms']))
        dp = max(abs(round(a, 3) - b['charge_penalty']) for a, b in zip(pen, mol['atoms']))
        ok = dq < 1e-9 and dp < 1e-9 and not missing
        ok_all &= ok
        print(f"{mol['name']}: max |dq| = {dq:.3f} e, max |d penalty| = {dp:.3f}, missing = {len(missing)} -> {'EXACT' if ok else 'MISMATCH'}")
    return 0 if ok_all else 1


def cmd_training(args):
    """Formula check on training residues: formal charges q0 are inferred as round(q_ref - q_increments)."""
    table = load_increments(f"{args.toppar}/par_all36_cgenff_increments.prm")
    res, cur = [], None
    for line in open(f"{args.toppar}/top_all36_cgenff.rtf", errors='ignore'):
        w = line.split('!')[0].split()
        if not w: continue
        if w[0] in ('RESI', 'PRES'):
            cur = {'net': float(w[2]), 'atoms': [], 'types': {}, 'q': {}, 'bonds': [], 'skip': w[0] == 'PRES'}; res.append(cur)
        elif cur is None: continue
        elif w[0] == 'ATOM': cur['atoms'].append(w[1]); cur['types'][w[1]] = w[2]; cur['q'][w[1]] = float(w[3])
        elif w[0] in ('BOND', 'DOUBLE', 'TRIPLE'):
            for a, b in zip(w[1::2], w[2::2]):
                if a[0] in '+-' or b[0] in '+-': cur['skip'] = True
                cur['bonds'].append((a, b))
        elif w[0] in ('LONEPAIR', 'LONE'): cur['skip'] = True
        elif w[0] == 'END': break
    errs, per_mol, used = [], [], 0
    for r in res:
        if r['skip'] or len(r['atoms']) < 2: continue
        idx = {a: i for i, a in enumerate(r['atoms'])}
        bonds = sorted({tuple(sorted((idx[a], idx[b]))) for a, b in r['bonds'] if a in idx and b in idx})
        types = [r['types'][a] for a in r['atoms']]
        q, _, missing = assign(types, bonds, [0] * len(types), table)
        if missing or any(_lookup(table, tuple(types[a] for a in inst))[0] is None for inst in _instances(len(types), bonds)):
            continue
        used += 1
        e = [r['q'][a] - q[i] - round(r['q'][a] - q[i]) for i, a in enumerate(r['atoms'])]
        errs += e; per_mol.append(math.sqrt(sum(x * x for x in e) / len(e)))
    print(f"training residues with complete increments: {used}")
    print(f"median per-molecule RMSE {statistics.median(per_mol):.4f} e; atoms within 0.01 e: "
          f"{sum(abs(x) < 0.01 for x in errs) / len(errs):.1%}; within 0.05 e: {sum(abs(x) < 0.05 for x in errs) / len(errs):.1%}")
    print("note: residual is the fit error of the increment model vs hand-optimized topology charges; the program's "
          "output for these molecules is the increment value, not the topology charge.")
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fixtures'); f.add_argument('--toppar', required=True); f.add_argument('--fixtures', required=True)
    t = sub.add_parser('training'); t.add_argument('--toppar', required=True)
    a = ap.parse_args()
    sys.exit({'fixtures': cmd_fixtures, 'training': cmd_training}[a.cmd](a))

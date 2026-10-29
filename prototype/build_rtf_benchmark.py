"""Benchmark B1: typing benchmark from the public CGenFF 5.0 topology file (top_all36_cgenff.rtf).
For each residue with a SMILES comment, map RTF atoms 1:1 onto an RDKit molecule (explicit H, bond orders)
by graph isomorphism, keep it only if the bond graphs are identical, and store the reference atom types.
Residues whose SMILES comment has lost bond orders (radicals) are repaired from connectivity + net charge.
Output: <out>.pkl = list of dicts {name, origin, smiles, mol (binary), types (RDKit atom order), split}
        <out>.split.json = {train: [...], test: [...]} (fixed seed 42; 80/20)
usage: python build_rtf_benchmark.py <top_all36_cgenff.rtf> <out.pkl>"""
import sys, pickle, random, json
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

rtf, out = sys.argv[1], sys.argv[2]
ELEM = [('HG', 'H'), ('CG', 'C'), ('NG', 'N'), ('OG', 'O'), ('SG', 'S'), ('PG', 'P'), ('FG', 'F'), ('CLG', 'Cl'),
        ('BRG', 'Br'), ('IG', 'I'), ('BG', 'B'), ('SEG', 'Se'), ('ALG', 'Al')]

def elem(t):
    for k, v in ELEM:
        if t.startswith(k): return v
    return None

res, cur = [], None
for l in open(rtf, errors='ignore'):
    s = l.split('!')[0].split()
    if l.startswith('RESI'):
        cur = {'name': l.split()[1], 'net': float(l.split()[2]), 'atoms': [], 'types': {}, 'bonds': [], 'smiles': None, 'lp': False}
        res.append(cur); continue
    if l.startswith('PRES'): cur = None; continue
    if cur is None: continue
    if 'SMILES:' in l and cur['smiles'] is None: cur['smiles'] = l.split('SMILES:')[1].split()[0]
    if not s: continue
    if s[0] == 'ATOM': cur['atoms'].append(s[1]); cur['types'][s[1]] = s[2]
    elif s[0] in ('BOND', 'DOUBLE', 'TRIPLE'): cur['bonds'] += list(zip(s[1::2], s[2::2]))
    elif s[0] in ('LONEPAIR', 'LONE'): cur['lp'] = True
    elif s[0] == 'END': break

from rdkit.Chem import rdDetermineBonds
def repair(mol, net):
    m = Chem.RWMol(mol)
    for a in m.GetAtoms(): a.SetFormalCharge(0); a.SetNumRadicalElectrons(0); a.SetNoImplicit(True); a.SetIsAromatic(False)
    for b in m.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    try:
        rdDetermineBonds.DetermineBondOrders(m, charge=int(round(net)), embedChiral=False)
        m = m.GetMol(); Chem.SanitizeMol(m)
        if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()): return None
        return m
    except Exception:
        return None

nrad = [0, 0]
data, skipped = [], 0
for r in res:
    if not r['smiles']: continue
    mol = Chem.MolFromSmiles(r['smiles'])
    if mol is None: skipped += 1; continue
    mol = Chem.AddHs(mol)
    names = [a for a in r['atoms'] if elem(r['types'][a])]
    if len(names) != mol.GetNumAtoms(): skipped += 1; continue
    idx = {n: i for i, n in enumerate(names)}
    q = Chem.RWMol()
    for n in names: q.AddAtom(Chem.AtomFromSmarts(f"[#{Chem.Atom(elem(r['types'][n])).GetAtomicNum()}]"))
    bonds = set()
    for a, b in r['bonds']:
        if a in idx and b in idx and frozenset((a, b)) not in bonds:
            bonds.add(frozenset((a, b)))
            q.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
            q.ReplaceBond(q.GetBondBetweenAtoms(idx[a], idx[b]).GetIdx(), Chem.BondFromSmarts('~'))
    if len(bonds) != mol.GetNumBonds(): skipped += 1; continue
    mt = mol.GetSubstructMatch(q.GetMol())
    if len(mt) != mol.GetNumAtoms(): skipped += 1; continue
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        fixed = repair(mol, r['net'])
        if fixed is None: skipped += 1; nrad[1] += 1; continue
        mol = fixed; nrad[0] += 1
    types = [None] * mol.GetNumAtoms()
    for qi, mi in enumerate(mt): types[mi] = r['types'][names[qi]]
    data.append({'name': r['name'], 'origin': 'zinc' if r['name'].startswith('230') else 'legacy',
                 'smiles': r['smiles'], 'mol': mol.ToBinary(), 'types': types, 'lp': r['lp']})

random.seed(42)
random.shuffle(data)
n_test = len(data) // 5
for i, d in enumerate(data): d['split'] = 'test' if i < n_test else 'train'
pickle.dump(data, open(out, 'wb'))
json.dump({s: sorted(d['name'] for d in data if d['split'] == s) for s in ('train', 'test')}, open(out.rsplit('.', 1)[0] + '.split.json', 'w'), indent=0)
print(f"radical SMILES repaired: {nrad[0]}, unrepairable: {nrad[1]}")
print(f"molecules: {len(data)} (skipped {skipped}); train {len(data)-n_test}, held-out {n_test}; "
      f"atoms {sum(len(d['types']) for d in data)}; zinc {sum(d['origin']=='zinc' for d in data)}, legacy {sum(d['origin']=='legacy' for d in data)}")

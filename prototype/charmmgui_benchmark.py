"""Benchmark B3: CHARMM-GUI GROMACS topologies ([A].itp) with reference CGenFF types (CGenFF <= 4.6 era).
Builds RDKit molecules (elements from masses, bond orders inferred from connectivity + net charge),
then scores a typer module and reports per-molecule results.
usage: python charmmgui_benchmark.py <itp_dir> <typer.py> [--detail] [--out ext_dataset.pkl]"""
import sys, os, glob, importlib.util, collections, pickle
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds
RDLogger.DisableLog('rdApp.*')

itp_dir, typer_path = sys.argv[1], sys.argv[2]
detail = '--detail' in sys.argv
spec = importlib.util.spec_from_file_location('typer', typer_path); typer = importlib.util.module_from_spec(spec); spec.loader.exec_module(typer)
MASS = [(1.008, 'H'), (10.811, 'B'), (12.011, 'C'), (14.007, 'N'), (15.9994, 'O'), (18.998, 'F'), (30.974, 'P'), (32.06, 'S'),
        (35.45, 'Cl'), (79.904, 'Br'), (126.904, 'I'), (196.97, 'Au'), (22.99, 'Na'), (78.96, 'Se')]
ALT = {'CG2DC1': 'CG2DC?', 'CG2DC2': 'CG2DC?', 'CG2D1O': 'CG2D?O', 'CG2D2O': 'CG2D?O',
       'CG25C1': 'CG25C?', 'CG25C2': 'CG25C?', 'CG251O': 'CG25?O', 'CG252O': 'CG25?O'}
canon = lambda t: ALT.get(t, t)


def elem(mass):
    return min(MASS, key=lambda x: abs(x[0] - mass))[1] if mass > 0.5 else None


def parse(fn):
    sec = None; atoms = []; bonds = []; name = None
    for l in open(fn, errors='ignore'):
        s = l.split(';')[0].strip()
        if not s: continue
        if s.startswith('['): sec = s.strip('[] ').lower(); continue
        w = s.split()
        if sec == 'moleculetype' and name is None: name = w[0]
        elif sec == 'atoms' and len(w) >= 8: atoms.append({'type': w[1], 'name': w[4], 'q': float(w[6]), 'mass': float(w[7])})
        elif sec == 'bonds' and len(w) >= 2: bonds.append((int(w[0]) - 1, int(w[1]) - 1))
    return name, atoms, bonds


def build(atoms, bonds):
    keep = [i for i, a in enumerate(atoms) if elem(a['mass'])]          # drop virtual sites / lone pairs (mass 0)
    idx = {i: k for k, i in enumerate(keep)}
    m = Chem.RWMol()
    for i in keep:
        a = Chem.Atom(elem(atoms[i]['mass'])); a.SetNoImplicit(True); m.AddAtom(a)
    for i, j in bonds:
        if i in idx and j in idx and m.GetBondBetweenAtoms(idx[i], idx[j]) is None:
            m.AddBond(idx[i], idx[j], Chem.BondType.SINGLE)
    net = round(sum(atoms[i]['q'] for i in range(len(atoms))))
    if net == 0 and any(atoms[i]['type'] == 'SG302' for i in keep): net = -1   # thiolate with custom neutral charges (MBN)
    rdDetermineBonds.DetermineBondOrders(m, charge=int(net), embedChiral=False)
    mol = m.GetMol(); Chem.SanitizeMol(mol)
    return mol, keep, net


rows = []; data = []
for fn in sorted(glob.glob(os.path.join(itp_dir, '*.itp'))):
    base = os.path.basename(fn)
    if base.startswith('cgen_'): continue
    name, atoms, bonds = parse(fn)
    if not atoms and os.path.exists(os.path.join(itp_dir, 'cgen_' + base)):   # files swapped (MBN)
        name, atoms, bonds = parse(os.path.join(itp_dir, 'cgen_' + base))
    ref = [a['type'] for a in atoms]
    if not atoms:
        rows.append((base, name, 0, None, None, 'no atoms section', ref)); continue
    try:
        mol, keep, net = build(atoms, bonds)
    except Exception as e:
        rows.append((base, name, len(atoms), None, None, f'could not build molecule ({type(e).__name__})', ref)); continue
    if mol.GetNumAtoms() < 2 or (all(a.GetSymbol() in ('H', 'O') for a in mol.GetAtoms()) and mol.GetNumAtoms() <= 3):
        rows.append((base, name, len(atoms), None, None, 'water / single atom (not a CGenFF ligand)', ref)); continue
    reft = [ref[i] for i in keep]
    if any(t[:2] in ('CC', 'HC', 'OC') and not t.startswith(('CG', 'HG', 'OG', 'CL')) for t in reft):
        rows.append((base, name, len(keep), None, None, 'CHARMM36 carbohydrate types, not CGenFF', reft)); continue
    data.append({'name': name, 'file': base, 'mol': mol.ToBinary(), 'types': reft, 'charges': [atoms[i]['q'] for i in keep], 'net': net,
                 'smiles': Chem.MolToSmiles(Chem.RemoveHs(mol))})
    try:
        pred = typer.type_molecule(mol)
    except Exception as e:
        rows.append((base, name, len(keep), None, None, f'typer crashed: {e}', reft)); continue
    bad = [(mol.GetAtomWithIdx(i).GetSymbol() + str(i), t, p) for i, (t, p) in enumerate(zip(reft, pred)) if canon(t) != canon(p)]
    rows.append((base, name, len(keep), len(keep) - len(bad), bad, Chem.MolToSmiles(Chem.RemoveHs(mol)), reft))

if '--out' in sys.argv: pickle.dump(data, open(sys.argv[sys.argv.index('--out') + 1], 'wb'))
scored = [r for r in rows if r[3] is not None]
na = sum(r[2] for r in scored); ok = sum(r[3] for r in scored)
print(f"species: {len(rows)}; scored: {len(scored)}; atoms correct {ok}/{na} ({ok/na:.1%}); "
      f"molecules fully correct {sum(r[3]==r[2] for r in scored)}/{len(scored)}")
for r in rows:
    if r[3] is None: print(f"  SKIP {r[0]:14s} {r[5]}"); continue
    status = 'OK ' if r[3] == r[2] else 'ERR'
    print(f"  {status} {r[0]:14s} {r[3]:3d}/{r[2]:<3d} {r[5][:70]}")
    if detail and r[4]:
        for a, t, p in r[4][:8]: print(f"        {a:6s} ref {t:7s} pred {p}")

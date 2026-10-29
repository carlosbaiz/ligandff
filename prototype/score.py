"""Score a typer module against the benchmark.
usage: python score.py dataset.pkl typer.py [--examples N] [--split train|test|both]
The typer module must define type_molecule(mol) -> list of type strings (RDKit atom order).
Conjugated-chain types that differ only in the arbitrary 1/2 'altnum' digit are compared as equal."""
import sys, pickle, collections, importlib.util, re, traceback
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

ds, typer_path = sys.argv[1], sys.argv[2]
nex = int(sys.argv[sys.argv.index('--examples') + 1]) if '--examples' in sys.argv else 0
show = sys.argv[sys.argv.index('--split') + 1] if '--split' in sys.argv else 'train'
spec = importlib.util.spec_from_file_location('typer', typer_path); typer = importlib.util.module_from_spec(spec); spec.loader.exec_module(typer)
data = pickle.load(open(ds, 'rb'))

ALT = {'CG2DC1': 'CG2DC?', 'CG2DC2': 'CG2DC?', 'CG2D1O': 'CG2D?O', 'CG2D2O': 'CG2D?O',
       'CG25C1': 'CG25C?', 'CG25C2': 'CG25C?', 'CG251O': 'CG25?O', 'CG252O': 'CG25?O'}
canon = lambda t: ALT.get(t, t)

stats = {s: {'atoms': 0, 'ok': 0, 'mols': 0, 'molok': 0, 'crash': 0} for s in ('train', 'test')}
conf = collections.Counter(); ex = collections.defaultdict(list); byorigin = collections.defaultdict(lambda: [0, 0])
for d in data:
    mol = Chem.Mol(d['mol']); st = stats[d['split']]; st['mols'] += 1
    try:
        pred = typer.type_molecule(mol)
    except Exception as e:
        st['crash'] += 1; st['atoms'] += len(d['types'])
        if d['split'] == show and len(ex['CRASH']) < 3: ex['CRASH'].append((d['name'], d['smiles'], traceback.format_exc(limit=2)))
        continue
    good = True
    for i, (t, p) in enumerate(zip(d['types'], pred)):
        st['atoms'] += 1
        if canon(t) == canon(p): st['ok'] += 1
        else:
            good = False
            if d['split'] == show or show == 'both':
                conf[(t, p)] += 1
                if len(ex[(t, p)]) < 2:
                    m2 = Chem.Mol(mol); m2.GetAtomWithIdx(i).SetAtomMapNum(99)
                    ex[(t, p)].append((d['name'], Chem.MolToSmiles(Chem.RemoveHs(m2, sanitize=False)) if t[0] != 'H' else Chem.MolToSmiles(m2)))
    st['molok'] += good
    if d['split'] == 'test': byorigin[d['origin']][0] += good; byorigin[d['origin']][1] += 1

for s in ('train', 'test'):
    x = stats[s]
    print(f"{s:5s}: atoms {x['ok']/x['atoms']:.2%}  molecules fully correct {x['molok']/x['mols']:.1%} "
          f"({x['molok']}/{x['mols']})  crashes {x['crash']}")
print("held-out molecules fully correct by origin:", {k: f"{v[0]/v[1]:.1%} of {v[1]}" for k, v in byorigin.items()})
if nex:
    print(f"\ntop confusions on {show} (true -> predicted), with example SMILES (atom of interest marked :99):")
    for (t, p), n in conf.most_common(nex):
        print(f"  {n:4d}  {t:7s} -> {p}   e.g. {ex[(t, p)][0][1]}  [{ex[(t, p)][0][0]}]")
    for c in ex['CRASH']: print('CRASH', c[0], c[1], c[2][-400:])

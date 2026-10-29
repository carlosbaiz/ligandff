"""Query the TRAINING split: for atoms of given true type(s), tabulate a feature of their environment.
usage: python query.py dataset.pkl TYPE[,TYPE] feature   (feature: parent | nbrs | nbrtypes | ring)"""
import sys, pickle, collections
from rdkit import Chem
data = pickle.load(open(sys.argv[1], 'rb')); want = set(sys.argv[2].split(',')); feat = sys.argv[3]
c = collections.Counter(); ex = {}
for d in data:
    if d['split'] != 'train': continue
    mol = Chem.Mol(d['mol']); T = d['types']
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if T[i] not in want: continue
        nb = list(a.GetNeighbors())
        if feat == 'parent': k = (T[i], T[nb[0].GetIdx()])
        elif feat == 'nbrtypes': k = (T[i], tuple(sorted(T[n.GetIdx()] for n in nb)))
        elif feat == 'nbrs': k = (T[i], a.GetTotalNumHs(includeNeighbors=True), tuple(sorted(f"{n.GetSymbol()}{'ar' if n.GetIsAromatic() else ''}{n.GetFormalCharge() or ''}" for n in nb if n.GetSymbol() != 'H')))
        elif feat == 'ring': k = (T[i], tuple(sorted(len(r) for r in mol.GetRingInfo().AtomRings() if i in r)), a.GetIsAromatic())
        c[k] += 1
        if k not in ex:
            m2 = Chem.Mol(mol); m2.GetAtomWithIdx(i).SetAtomMapNum(99); ex[k] = Chem.MolToSmiles(Chem.RemoveHs(m2, sanitize=False)) if T[i][0] != 'H' else d['smiles']
for k, n in c.most_common(int(sys.argv[4]) if len(sys.argv) > 4 else 25):
    print(f"{n:5d}  {k}   e.g. {ex[k][:90]}")

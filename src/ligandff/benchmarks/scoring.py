"""Typing metrics (SPEC 7.3). Per-molecule details are produced only for development splits (SPEC 7.2)."""
import collections
import time
import traceback

from rdkit import Chem

from .common import canonical_type

DEV_SPLITS = ("train", "dev")
MAX_EXAMPLES = 2


def _marked_smiles(mol, atom_index, is_h):
    m = Chem.Mol(mol)
    m.GetAtomWithIdx(atom_index).SetAtomMapNum(99)
    return Chem.MolToSmiles(m) if is_h else Chem.MolToSmiles(Chem.RemoveHs(m, sanitize=False))


def score_typing(cases, type_molecule, n_confusions=30):
    """Score a typer on cases; returns {split: metrics}. Confusions and failures only for DEV_SPLITS."""
    by_split = collections.defaultdict(list)
    for c in cases:
        by_split[c.split].append(c)
    return {split: _score_split(rows, type_molecule, split in DEV_SPLITS, n_confusions)
            for split, rows in sorted(by_split.items())}


def _score_split(cases, type_molecule, details, n_confusions):
    n_atoms = n_ok_atoms = n_ok_mols = n_crash = 0
    tp, fp, fn = collections.Counter(), collections.Counter(), collections.Counter()
    confusions, examples, failures = collections.Counter(), collections.defaultdict(list), []
    t0 = time.perf_counter()
    for c in cases:
        n_atoms += len(c.types)
        try:
            pred = type_molecule(Chem.Mol(c.mol))
        except Exception:
            n_crash += 1
            fn.update(canonical_type(t) for t in c.types)
            if details:
                failures.append({"name": c.name, "smiles": c.smiles, "crash": traceback.format_exc(limit=2)[-400:]})
            continue
        wrong = []
        for i, (t, p) in enumerate(zip(c.types, pred)):
            t, p = canonical_type(t), canonical_type(p)
            if t == p:
                n_ok_atoms += 1; tp[t] += 1
                continue
            fn[t] += 1; fp[p] += 1
            wrong.append((i, t, p))
            if details:
                confusions[(t, p)] += 1
                if len(examples[(t, p)]) < MAX_EXAMPLES:
                    examples[(t, p)].append(
                        {"name": c.name, "smiles": _marked_smiles(c.mol, i, t.startswith("H"))})
        n_ok_mols += not wrong
        if details and wrong:
            failures.append({"name": c.name, "smiles": c.smiles,
                             "atoms": [{"index": i, "ref": t, "pred": p} for i, t, p in wrong]})
    elapsed = time.perf_counter() - t0
    types = sorted(set(tp) | set(fp) | set(fn))
    out = {
        "molecules": len(cases), "molecules_correct": n_ok_mols, "molecule_accuracy": n_ok_mols / len(cases),
        "atoms": n_atoms, "atoms_correct": n_ok_atoms, "atom_accuracy": n_ok_atoms / n_atoms, "crashes": n_crash,
        "molecules_per_second": len(cases) / elapsed if elapsed else None,
        "per_type": {t: {"tp": tp[t], "fp": fp[t], "fn": fn[t],
                         "precision": tp[t] / (tp[t] + fp[t]) if tp[t] + fp[t] else None,
                         "recall": tp[t] / (tp[t] + fn[t]) if tp[t] + fn[t] else None} for t in types},
    }
    if details:
        out["confusions"] = [{"ref": t, "pred": p, "count": n, "examples": examples[(t, p)]}
                             for (t, p), n in confusions.most_common(n_confusions)]
        out["failures"] = failures
    return out

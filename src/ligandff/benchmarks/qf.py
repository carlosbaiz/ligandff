"""Set QF: typing-rule formal charges vs the training charges (SPEC 5.7).

For every B1/B2 molecule without lone pairs whose increments are complete and whose atom types the typer gets right,
the rule charge of each atom must equal q0 = round(q_ref - q_increments), the integer formal charge implied by the
hand-optimized topology charges (FINDINGS 2.2). Per-atom mismatches are listed for development splits only.
"""
import collections

from ..charges import assign_charges, parameter_instances
from .common import canonical_type
from .scoring import DEV_SPLITS


def inferred_formal_charges(case, table):
    """q0 per atom from the reference types and charges, or None if increments are incomplete."""
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in case.mol.GetBonds()]
    result = assign_charges(case.types, bonds, [0] * len(case.types), table)
    if result.missing or any(table.lookup([case.types[a] for a in inst])[0] is None
                             for inst in parameter_instances(len(case.types), bonds)):
        return None
    return [round(q_ref - q) for q_ref, q in zip(case.meta["charges"], result.charges)]


def score_formal_charges(sets, type_molecule, table, n_examples=30):
    """sets: {name: cases}; type_molecule(rdkit_mol) -> TypedMolecule (types and formal_charges)."""
    out = {}
    for name, cases in sets.items():
        per_split = collections.defaultdict(lambda: {"molecules": 0, "correct": 0, "atoms_wrong": 0})
        wrong = collections.Counter()
        examples = {}
        for case in cases:
            if case.meta.get("lonepair") or case.meta.get("charges") is None:
                continue
            q0 = inferred_formal_charges(case, table)
            if q0 is None:
                continue
            typed = type_molecule(case.mol)
            if any(canonical_type(a) != canonical_type(b) for a, b in zip(typed.types, case.types)):
                continue
            s = per_split[case.split]
            s["molecules"] += 1
            bad = [i for i, (a, b) in enumerate(zip(typed.formal_charges, q0)) if a != b]
            s["correct"] += not bad
            s["atoms_wrong"] += len(bad)
            if case.split in DEV_SPLITS:
                for i in bad:
                    key = (case.types[i], q0[i], typed.formal_charges[i])
                    wrong[key] += 1
                    examples.setdefault(key, case.name)
        out[name] = {"splits": dict(per_split),
                     "mismatches": [{"type": t, "inferred": q, "rules": r, "count": n, "example": examples[(t, q, r)]}
                                    for (t, q, r), n in wrong.most_common(n_examples)]}
    return out

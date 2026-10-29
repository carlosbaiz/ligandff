"""Formula check on the CGenFF training residues (FINDINGS 2.2).

For every RESI of top_all36_cgenff.rtf without lone pairs or inter-residue bonds whose parameters all have
increments (no N/A), increment charges + an inferred integer formal charge q0 = round(q_ref - q_increments) are
compared with the hand-optimized topology charges. The residual is the fit error of the increment model, not an
error of the engine: for these molecules the program outputs the increment charges.
"""
import collections
import math
import statistics
from dataclasses import dataclass, field

from ..io.rtf import read_rtf
from ..paths import TOPOLOGY_FILE, toppar_dir
from .increments import IncrementTable, assign_charges, parameter_instances


@dataclass
class TrainingCheck:
    residues: int
    median_rmse: float
    fraction_within_001: float
    fraction_within_005: float
    formal_charges: dict = field(default_factory=dict)    # (type, inferred q0) -> count, q0 != 0 only

    def as_dict(self):
        return {"residues": self.residues, "median_rmse": self.median_rmse,
                "fraction_within_0.01": self.fraction_within_001, "fraction_within_0.05": self.fraction_within_005,
                "formal_charges": [{"type": t, "q0": q, "count": n} for (t, q), n in
                                   sorted(self.formal_charges.items(), key=lambda x: -x[1])]}


def _usable(res):
    return res.kind == "RESI" and not res.has_lonepair and not res.has_external_bond and len(res.atoms) >= 2


def training_check(toppar=None, table=None):
    toppar = toppar_dir(toppar)
    table = table or IncrementTable.from_toppar(toppar)
    errors, rmse, q0_counts = [], [], collections.Counter()
    for res in read_rtf(toppar / TOPOLOGY_FILE):
        if not _usable(res):
            continue
        idx = {a.name: i for i, a in enumerate(res.atoms)}
        bonds = sorted({tuple(sorted((idx[a], idx[b]))) for a, b in res.bonds if a in idx and b in idx})
        types = [a.type for a in res.atoms]
        result = assign_charges(types, bonds, [0] * len(types), table)
        if result.missing or any(table.lookup([types[a] for a in inst])[0] is None
                                 for inst in parameter_instances(len(types), bonds)):
            continue
        err = []
        for atom, q in zip(res.atoms, result.charges):
            q0 = round(atom.charge - q)
            err.append(atom.charge - q - q0)
            if q0:
                q0_counts[(atom.type, q0)] += 1
        errors += err
        rmse.append(math.sqrt(sum(x * x for x in err) / len(err)))
    return TrainingCheck(len(rmse), statistics.median(rmse), sum(abs(x) < 0.01 for x in errors) / len(errors),
                         sum(abs(x) < 0.05 for x in errors) / len(errors), dict(q0_counts))

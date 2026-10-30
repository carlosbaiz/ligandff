"""OUT: end-to-end output check (SPEC section 1, 5.14): every B1/B2/B3 molecule goes through the full pipeline
(typing, parameters by analogy, charges, lone pairs), is written as `.str` + PSF and read back by OpenMM's CHARMM
readers together with the CGenFF files; the system must have one particle per atom and lone pair and a finite energy.

Not in the default sets (needs OpenMM; ~2 min on 14 cores). Reports pipeline throughput (SPEC 7.3) and the penalty
distribution. Per-molecule failures are listed for development splits only (SPEC 7.2).
"""
import collections
import math
import time

from ..molecule import InputMolecule
from ..params import default_engine
from ..pipeline import increment_table, parameterize

RESNAME = "MOL"
PENALTY_BINS = ((10, "< 10"), (50, "10-50"), (math.inf, "> 50"))     # SPEC 5.10 guidance on the RESI penalty
DEV_SPLITS = ("train", "dev")
_TOPPAR = None


def _bin(x):
    return next(label for limit, label in PENALTY_BINS if x < limit)


def check_molecule(case, toppar=None):
    """One molecule through pipeline + OpenMM; a dict of results (never raises)."""
    from ..export import openmm_system, potential_energy
    row = {"name": case.name, "errors": [], "built": False, "finite": False}
    t = time.perf_counter()
    try:
        r = parameterize(InputMolecule.from_rdkit(case.mol, case.name), RESNAME, toppar)
    except Exception as exc:                                         # pipeline crash
        row["errors"].append(f"pipeline: {exc!r}"[:200])
        return row
    row["pipeline_seconds"] = time.perf_counter() - t
    row.update(lone_pairs=len(r.lone_pairs), analogs=len(r.analogs), param_penalty=r.param_penalty,
               charge_penalty=r.charge_penalty, missing=len(r.assignment.missing))
    row["errors"] += r.errors
    try:
        system, _, pos = openmm_system(r, toppar)
        row["built"] = system.getNumParticles() == len(r.atoms) + len(r.lone_pairs)
        row["finite"] = math.isfinite(potential_energy(system, pos))
    except Exception as exc:
        row["errors"].append(f"openmm: {exc!r}"[:200])
    return row


def _check_chunk(chunk):
    return [check_molecule(case, _TOPPAR) for case in chunk]


def _rows(cases, toppar, processes):
    global _TOPPAR
    import multiprocessing as mp
    _TOPPAR = toppar
    default_engine(toppar), increment_table(toppar)                  # load once, shared by the forked workers
    if processes <= 1 or "fork" not in mp.get_all_start_methods():
        return _check_chunk(cases)
    chunks = [cases[i::processes] for i in range(processes)]
    with mp.get_context("fork").Pool(processes) as pool:
        parts = pool.map(_check_chunk, chunks)
    out = [None] * len(cases)
    for i, part in enumerate(parts):                                 # chunk i holds cases i, i + processes, ...
        for k, row in enumerate(part):
            out[i + k * processes] = row
    return out


def score_output(cases_by_set, toppar=None, processes=None):
    """{set: {split: counts}} plus failures (development splits) and pipeline throughput (one core)."""
    import os
    processes = processes or min(os.cpu_count() or 1, 16)
    flat = [(name, case) for name, cases in cases_by_set.items() for case in cases]
    rows = _rows([case for _, case in flat], toppar, processes)
    per = collections.defaultdict(collections.Counter)
    failures, seconds = [], []
    for (name, case), row in zip(flat, rows):
        s = per[(name, case.split)]
        ok = not row["errors"] and row["built"] and row["finite"]
        s["molecules"] += 1
        s["ok"] += ok
        s["pipeline_errors"] += "pipeline_seconds" not in row or any(not e.startswith("openmm") for e in row["errors"])
        s["openmm_built"] += row["built"] and row["finite"]
        if "pipeline_seconds" in row:
            seconds.append(row["pipeline_seconds"])
            s["with_lone_pairs"] += row["lone_pairs"] > 0
            s["with_analogs"] += row["analogs"] > 0
            s["missing_parameters"] += row["missing"] > 0
            s[f"param_penalty {_bin(row['param_penalty'])}"] += 1
            s[f"charge_penalty {_bin(row['charge_penalty'])}"] += 1
        if not ok and case.split in DEV_SPLITS:
            failures.append({"set": name, "name": case.name, "errors": row["errors"]})
    splits = {f"{n}/{sp}": dict(c) for (n, sp), c in sorted(per.items())}
    total = sum(seconds)
    return {"splits": splits, "failures": failures,
            "molecules_per_second": len(seconds) / total if total else 0.0, "processes": processes}

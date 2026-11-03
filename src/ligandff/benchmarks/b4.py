"""B4: the program's own analogy choices recorded in the parameter file (SPEC 7.1, 5.10).

Each bonded parameter annotated `from <source>, penalty= <tps>` is a case: the missing types, the source the program
chose and its TPS. Leave-one-out (SPEC 5.10 step 3): search an analog for the case's types among the other parameters
and compare source and TPS.

Candidate pool: 6,952 of the 7,234 annotations come from one batch (comment starting `aa,`). None of them uses a
parameter with an `aa,` comment as its source, annotated or fitted, although those are 18,852 of the 28,199 parameters
(67 %, FINDINGS §10): the batch was generated against a database without any `aa,` parameter. For a batch case the
pool therefore excludes every `aa,` parameter; for the other cases it excludes only the case itself. This replays the
conditions of the recorded choice; the shipped tool searches the full file. `score` also reports agreement by origin
(batch / non-batch) for this pool and for the full file (every other parameter a candidate), as aggregates.

Split 80/20 with seed 42 (`benchmarks/b4_annotations.split.json`); fitting may look at B4-train only, B4-test is
reported as aggregates (SPEC 7.2).
"""
import collections
import json
from dataclasses import dataclass

from ..paths import REPO_ROOT

MANIFEST = REPO_ROOT / "benchmarks" / "b4_annotations.split.json"
BATCH_PREFIX = "aa,"


@dataclass(frozen=True)
class AnalogyCase:
    kind: str
    types: tuple
    source: tuple            # as written in the annotation
    penalty: float
    batch: bool              # member of the one-run batch (comment starts with "aa,")
    split: str = "unassigned"

    @property
    def name(self):
        return f"{self.kind}:{' '.join(self.types)}"


def in_batch(prm):
    """Annotated member of the batch."""
    return prm.annotation is not None and prm.comment.startswith(BATCH_PREFIX)


def batch_era(prm):
    """Any parameter of the batch's project (annotated or fitted): absent when the batch analogies were made."""
    return prm.comment.startswith(BATCH_PREFIX)


def build(db, manifest=MANIFEST):
    """All annotated parameters of a ParameterDatabase as AnalogyCase, with the manifest's split."""
    split_of = {}
    if manifest is not None:
        for split, names in json.loads(open(manifest).read()).items():
            split_of.update({n: split for n in names})
    cases = []
    for kind, table in db.params.items():
        for prm in table.values():
            if prm.annotation is None:
                continue
            src, pen = prm.annotation
            c = AnalogyCase(kind, prm.types, tuple(src), pen, in_batch(prm))
            cases.append(AnalogyCase(c.kind, c.types, c.source, c.penalty, c.batch, split_of.get(c.name, "unassigned")))
    if manifest is not None and {c.name for c in cases} != set(split_of):
        raise ValueError(f"B4 build does not match manifest {manifest}")
    return cases


def _same_source(kind, found, source):
    if found in (source, source[::-1]):
        return True
    return kind == "improper" and found[0] == source[0] and sorted(found[1:]) == sorted(source[1:])


def leave_one_out(engine, case, full=False):
    """AnalogResult of the search for `case` with the B4 candidate pool (full=True: every other parameter), or None."""
    target = case.types

    def accept(prm):
        if prm.types == target:
            return False
        return full or not (case.batch and batch_era(prm))
    return engine.find(case.kind, target, accept=accept)


_ENGINE = None
_FULL = False


def _search(chunk):
    out = []
    for case in chunk:
        res = leave_one_out(_ENGINE, case, _FULL)
        out.append((case, res and (res.source, res.penalty)))
    return out


def _results(engine, cases, processes, full=False):
    """[(case, (source, penalty) or None)] in input order; forked worker processes where available."""
    global _ENGINE, _FULL
    import multiprocessing as mp
    _ENGINE, _FULL = engine, full
    if processes <= 1 or "fork" not in mp.get_all_start_methods():
        return _search(cases)
    chunks = [cases[i::processes] for i in range(processes)]
    with mp.get_context("fork").Pool(processes) as pool:
        parts = pool.map(_search, chunks)
    out = [None] * len(cases)
    for i, part in enumerate(parts):                      # chunk i holds cases i, i + processes, ...
        for k, (_, res) in enumerate(part):
            out[i + k * processes] = (cases[i + k * processes], res)
    return out


def score(engine, cases, n_examples=0, processes=None):
    """Per split: cases, same source, |dTPS| <= 0.5, both; per kind for development splits. Examples of failures
    (train split only) when n_examples > 0. processes: worker processes (default: CPU count, at most 16)."""
    import os
    processes = processes or min(os.cpu_count() or 1, 16)
    per = collections.defaultdict(lambda: collections.Counter())
    examples = []
    origin = collections.defaultdict(lambda: collections.Counter())
    for pool, full in (("B4 pool", False), ("full file", True)):
        for case, res in _results(engine, cases, processes, full):
            o = origin[(pool, case.split, "batch" if case.batch else "non-batch")]
            o["cases"] += 1
            o["same_source"] += res is not None and _same_source(case.kind, res[0], case.source)
            o["tps_within_0.5"] += res is not None and abs(res[1] - case.penalty) <= 0.5 + 1e-9
    for case, res in _results(engine, cases, processes):
        same = res is not None and _same_source(case.kind, res[0], case.source)
        close = res is not None and abs(res[1] - case.penalty) <= 0.5 + 1e-9
        for key in (case.split, f"{case.split}/{case.kind}"):
            s = per[key]
            s["cases"] += 1
            s["same_source"] += same
            s["tps_within_0.5"] += close
            s["both"] += same and close
        if not same and case.split == "train" and len(examples) < n_examples:
            examples.append({"case": case.name, "reference": [" ".join(case.source), case.penalty],
                             "found": res and [" ".join(res[0]), res[1]]})
    out = {}
    for key, s in per.items():
        if "/" in key and key.split("/")[0] != "train":
            continue
        out[key] = {k: s[k] for k in ("cases", "same_source", "tps_within_0.5", "both")}
        out[key].update({f"{k}_fraction": s[k] / s["cases"] for k in ("same_source", "tps_within_0.5", "both")})
    by_origin = [{"pool": pool, "split": split, "origin": o, **dict(c)} for (pool, split, o), c in sorted(origin.items())]
    return {"splits": out, "by_origin": by_origin, "examples": examples}

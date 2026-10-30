"""Benchmark B0: exact CGenFF program 5.0 output for three molecules (P-v5 SI Table S10; SPEC 7.1).

A charge engine is a callable engine(types, bonds, formal_charges, analogs) -> (charges, charge_penalties, missing),
where analogs = {missing parameter types: (source types, penalty)} as printed in the reference output.
"""
import json
import math

from ..paths import REPO_ROOT

FIXTURES = REPO_ROOT / "tests" / "fixtures" / "si_v5_reference.json"


def load(path=FIXTURES):
    return json.loads(open(path).read())


def engine_inputs(mol):
    """(types, bonds, formal_charges, analogs) for one fixture molecule."""
    names = [a["name"] for a in mol["atoms"]]
    idx = {n: i for i, n in enumerate(names)}
    return ([a["type"] for a in mol["atoms"]], [(idx[a], idx[b]) for a, b in mol["bonds"]],
            [mol["formal_charges"].get(n, 0) for n in names],
            {tuple(x["missing"]): (tuple(x["source"]), x["penalty"]) for x in mol["analog_parameters"]})


def score_charges(fixtures, engine):
    """Per-molecule charge and penalty agreement (outputs rounded to 0.001 as the program prints them)."""
    rows = []
    for mol in fixtures:
        q, pen, missing = engine(*engine_inputs(mol))
        dq = [round(x, 3) - a["charge"] for x, a in zip(q, mol["atoms"])]
        dp = [round(x, 3) - a["charge_penalty"] for x, a in zip(pen, mol["atoms"])]
        d_header = round(max(pen), 3) - mol["charge_penalty"]
        max_dq, max_dp = max(map(abs, dq)), max(map(abs, dp))
        rows.append({"name": mol["name"], "atoms": len(q), "analog_parameters": len(mol["analog_parameters"]),
                     "missing": [list(m) for m in missing],
                     "max_abs_dq": max_dq, "rms_dq": math.sqrt(sum(x * x for x in dq) / len(dq)),
                     "max_abs_dpenalty": max_dp, "d_header_charge_penalty": d_header,
                     "exact": not missing and max_dq < 1e-9 and max_dp < 1e-9 and abs(d_header) < 1e-9})
    return {"molecules": len(rows), "exact": sum(r["exact"] for r in rows), "per_molecule": rows}

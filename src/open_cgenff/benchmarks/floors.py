"""Regression floors (SPEC 7.4). CI fails if a run falls below any of them.

PILOT_BASELINE is the frozen pilot (prototype/typer_pilot.py); production must beat it. TYPING_FLOORS ratchet to the
production typer's results at the end of M3, so any regression fails CI; lower one only with an explanation logged
in docs/RULES_CHANGELOG.md. Floors are exact molecule counts."""
from fractions import Fraction

PILOT_BASELINE = {  # pilot on its own perception; B2 as rebuilt in M3 (charge-separation fix, b2.py)
    ("B1", "train"): Fraction(874, 904),
    ("B1", "test"): Fraction(212, 226),
    ("B2", "train"): Fraction(633, 728),
    ("B2", "test"): Fraction(161, 182),
    ("B3", "dev"): Fraction(1),
}

TYPING_FLOORS = {  # (set, split) -> minimum fraction of molecules with all atoms typed correctly (M3 production)
    ("B1", "train"): Fraction(895, 904),
    ("B1", "test"): Fraction(220, 226),
    ("B2", "train"): Fraction(720, 728),
    ("B2", "test"): Fraction(178, 182),
    ("B3", "dev"): Fraction(1),
}

FORMAL_CHARGE_FLOORS = {  # QF: (set, split) -> minimum number of molecules with all rule formal charges correct
    ("B1", "train"): 887,
    ("B1", "test"): 218,
    ("B2", "train"): 669,
    ("B2", "test"): 166,
}

IMPROPER_FLOORS = {  # QI: (set, split) -> minimum number of molecules with all improper centres correct (M4)
    ("B1", "train"): 884,
    ("B1", "test"): 219,
    ("B2", "train"): 705,
    ("B2", "test"): 175,
}

ANALOGY_FLOORS = {  # B4: (split) -> minimum (same source, TPS within 0.5) case counts (M5 fitted model)
    "train": (5383, 5358),     # of 5,787
    "test": (1333, 1338),      # of 1,447
}

# OUT: every B1/B2/B3 molecule's output builds an OpenMM system (SPEC section 1; M6)
OUTPUT_ALL_BUILD = True

# QT: the increment model's fit to the training charges (FINDINGS 2.2); a change means the engine or the data changed
CHARGE_TRAINING = {"residues": 2255, "median_rmse": 0.0098}


def check(report):
    """Return a list of floor violations (empty = pass)."""
    problems = []
    b0 = report.get("B0")
    if b0 is not None and b0["exact"] != b0["molecules"]:
        problems.append(f"B0: {b0['exact']}/{b0['molecules']} molecules exact (floor: all)")
    o = report.get("OUT")
    if o is not None and OUTPUT_ALL_BUILD:
        for key, m in o["splits"].items():
            if m["ok"] != m["molecules"]:
                problems.append(f"OUT {key}: {m['ok']}/{m['molecules']} molecules build an OpenMM system (floor: all)")
    qt = report.get("QT")
    if qt is not None and (qt["residues"] != CHARGE_TRAINING["residues"]
                           or round(qt["median_rmse"], 4) != CHARGE_TRAINING["median_rmse"]):
        problems.append(f"QT: {qt['residues']} residues, median RMSE {qt['median_rmse']:.4f} e "
                        f"(expected {CHARGE_TRAINING['residues']}, {CHARGE_TRAINING['median_rmse']})")
    qf = report.get("QF")
    if qf is not None:
        for (name, split), floor in FORMAL_CHARGE_FLOORS.items():
            m = qf.get(name, {}).get("splits", {}).get(split)
            if m is None or m["correct"] < floor:
                problems.append(f"QF-{name}-{split}: {m and m['correct']} molecules with correct formal charges "
                                f"(floor {floor})")
    b4 = report.get("B4")
    if b4 is not None:
        for split, (same, close) in ANALOGY_FLOORS.items():
            m = b4["splits"].get(split)
            if m is None or m["same_source"] < same or m["tps_within_0.5"] < close:
                problems.append(f"B4-{split}: {m and m['same_source']} same source, {m and m['tps_within_0.5']} TPS "
                                f"within 0.5 (floors {same}, {close})")
    qi = report.get("QI")
    if qi is not None:
        for (name, split), floor in IMPROPER_FLOORS.items():
            m = qi.get(name, {}).get("splits", {}).get(split)
            if m is None or m["correct"] < floor:
                problems.append(f"QI-{name}-{split}: {m and m['correct']} molecules with correct improper centres "
                                f"(floor {floor})")
    for (name, split), floor in TYPING_FLOORS.items():
        m = report.get(name, {}).get("typing", {}).get(split)
        if m is None:
            if name in report:
                problems.append(f"{name}-{split}: split missing from report")
            continue
        if Fraction(m["molecules_correct"], m["molecules"]) < floor:
            problems.append(f"{name}-{split}: {m['molecules_correct']}/{m['molecules']} molecules correct "
                            f"(floor {float(floor):.2%})")
    return problems

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

# QT: the increment model's fit to the training charges (FINDINGS 2.2); a change means the engine or the data changed
CHARGE_TRAINING = {"residues": 2255, "median_rmse": 0.0098}


def check(report):
    """Return a list of floor violations (empty = pass)."""
    problems = []
    b0 = report.get("B0")
    if b0 is not None and b0["exact"] != b0["molecules"]:
        problems.append(f"B0: {b0['exact']}/{b0['molecules']} molecules exact (floor: all)")
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

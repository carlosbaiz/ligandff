"""Regression floors (SPEC 7.4): the frozen pilot's results. CI fails if a run falls below any of them.
Typing floors are exact molecule counts: 874/904 = 96.68 % (printed as 96.7 % in the spec), 212/226 = 93.81 %."""
from fractions import Fraction

TYPING_FLOORS = {  # (set, split) -> minimum fraction of molecules with all atoms typed correctly
    ("B1", "train"): Fraction(874, 904),
    ("B1", "test"): Fraction(212, 226),
    ("B3", "dev"): Fraction(1),
}


def check(report):
    """Return a list of floor violations (empty = pass)."""
    problems = []
    b0 = report.get("B0")
    if b0 is not None and b0["exact"] != b0["molecules"]:
        problems.append(f"B0: {b0['exact']}/{b0['molecules']} molecules exact (floor: all)")
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

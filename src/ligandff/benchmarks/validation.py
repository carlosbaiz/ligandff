"""Final validation (M7): the definition of done of SPEC section 1, checked on a benchmark report.

usage: python -m ligandff.benchmarks.validation [--toppar DIR] [--out DIR]
Runs every set (B0, B1, B2, B3, B4, QT, QF, QI, OUT; OUT needs OpenMM), writes report.json/report.md and
validation.md (one row per target) to --out, and exits 1 if a target is missed. Test splits: aggregates only.
"""
import argparse
import json
import sys
from pathlib import Path

from ..paths import REPO_ROOT
from . import floors
from .runner import DEFAULT_SETS, run, to_markdown

SETS = DEFAULT_SETS + ("OUT",)

# SPEC section 1: (set, split) -> (minimum fraction of molecules fully correct, minimum atom accuracy)
TYPING_TARGETS = {("B1", "test"): (0.96, 0.995), ("B2", "test"): (0.95, 0.995), ("B3", "dev"): (1.0, 1.0)}
ANALOGY_TARGET = {"same_source": 0.85, "tps_within_0.5": 0.80}          # B4-test (M5)


def _row(target, measured, ok):
    return {"target": target, "measured": measured, "met": bool(ok)}


def targets(report):
    """[{target, measured, met}] for the section 1 definition of done."""
    rows = []
    for (name, split), (mol, atom) in TYPING_TARGETS.items():
        m = report[name]["typing"][split]
        rows.append(_row(f"types {name}-{split}: molecules >= {mol:.0%}, atoms >= {atom:.1%}",
                         f"{m['molecule_accuracy']:.2%} ({m['molecules_correct']}/{m['molecules']}), "
                         f"{m['atom_accuracy']:.2%}",
                         m["molecule_accuracy"] >= mol and m["atom_accuracy"] >= atom))
    b0 = report["B0"]
    rows.append(_row("charges B0: exact (charges and penalties to 0.001 e)", f"{b0['exact']}/{b0['molecules']}",
                     b0["exact"] == b0["molecules"]))
    qt = report["QT"]
    rows.append(_row(f"charges QT unchanged ({floors.CHARGE_TRAINING['residues']} residues, median RMSE "
                     f"{floors.CHARGE_TRAINING['median_rmse']} e)",
                     f"{qt['residues']} residues, {qt['median_rmse']:.4f} e",
                     not [v for v in floors.check({"QT": qt})]))
    for key, floor, label in (("QF", floors.FORMAL_CHARGE_FLOORS, "formal charges QF"),
                              ("QI", floors.IMPROPER_FLOORS, "impropers QI")):
        got = {(n, sp): m["correct"] for n, r in report[key].items() for sp, m in r["splits"].items()}
        rows.append(_row(f"{label} >= {'/'.join(str(v) for v in floor.values())} (SPEC 7.4)",
                         "/".join(str(got.get(k, 0)) for k in floor), all(got.get(k, 0) >= v for k, v in floor.items())))
    m = report["B4"]["splits"]["test"]
    rows.append(_row(f"analogy B4-test: same source >= {ANALOGY_TARGET['same_source']:.0%}, TPS within 0.5 >= "
                     f"{ANALOGY_TARGET['tps_within_0.5']:.0%}",
                     f"{m['same_source_fraction']:.1%}, {m['tps_within_0.5_fraction']:.1%} ({m['cases']} cases)",
                     all(m[f"{k}_fraction"] >= v for k, v in ANALOGY_TARGET.items())))
    o = report["OUT"]["splits"]
    n, ok = sum(x["molecules"] for x in o.values()), sum(x["ok"] for x in o.values())
    rows.append(_row("output: `.str` + PSF build an OpenMM system for every B1/B2/B3 molecule", f"{ok}/{n}", ok == n))
    return rows


def to_validation_markdown(report, rows):
    out = ["# Validation (M7): definition of done, SPEC section 1", "", f"{report['tool']} · {report['date']}", "",
           "| target | measured | met |", "|---|---|---|"]
    out += [f"| {r['target']} | {r['measured']} | {'yes' if r['met'] else 'NO'} |" for r in rows]
    out += ["", f"Regression floors (SPEC 7.4): {'all met' if not report['floor_violations'] else 'VIOLATED'}", ""]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--toppar", help="CHARMM toppar directory (default: data/ from scripts/fetch_data.sh)")
    ap.add_argument("--out", default=str(REPO_ROOT / "build" / "report"), help="output directory")
    a = ap.parse_args(argv)
    report = run(SETS, toppar=a.toppar)
    rows = targets(report)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=1))
    (out / "report.md").write_text(to_markdown(report))
    md = to_validation_markdown(report, rows)
    (out / "validation.md").write_text(md)
    print(md)
    return 0 if all(r["met"] for r in rows) and not report["floor_violations"] else 1


if __name__ == "__main__":
    sys.exit(main())

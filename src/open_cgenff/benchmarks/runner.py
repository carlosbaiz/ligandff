"""Benchmark runner: scores typing and charges on B0/B1/B3 and writes JSON + markdown reports (SPEC 7.3).
Set QT is the charge formula check on the training residues (FINDINGS 2.2).

usage: python -m open_cgenff.benchmarks [--sets B0,B1,B3,QT] [--typer FILE] [--charge-engine FILE]
                                        [--toppar DIR] [--out DIR] [--check-floors]
Hygiene (SPEC 7.2): confusion pairs and per-molecule failures are reported for development splits only.
"""
import argparse
import datetime
import importlib.util
import json
import sys
from pathlib import Path

from .. import banner
from ..charges import IncrementTable, assign_charges, training_check
from ..paths import REPO_ROOT, TOPOLOGY_FILE, INCREMENTS_FILE, toppar_dir
from . import b0, b1, b3, floors
from .common import load_typer
from .scoring import score_typing

DEFAULT_TYPER = REPO_ROOT / "prototype" / "typer_pilot.py"            # until the production typer exists (M3)
DEFAULT_SETS = ("B0", "B1", "B3", "QT")
PRODUCTION = "open_cgenff.charges"


def load_charge_engine(path, toppar):
    """Production engine, or an adapter for a module file with load_increments(prm) and
    assign(types, bonds, q0, table, analogs) -> (charges, penalties, missing), e.g. prototype/charges.py."""
    if path in (None, PRODUCTION):
        table = IncrementTable.from_toppar(toppar)

        def engine(types, bonds, q0, analogs):
            r = assign_charges(types, bonds, q0, table, analogs)
            return r.charges, r.penalties, r.missing
        return engine
    path = Path(path).resolve()
    spec = importlib.util.spec_from_file_location(f"charges_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    table = module.load_increments(str(toppar / INCREMENTS_FILE))
    return lambda types, bonds, q0, analogs: module.assign(types, bonds, q0, table, analogs)


def run(sets=DEFAULT_SETS, typer=DEFAULT_TYPER, charge_engine=PRODUCTION, toppar=None,
        n_confusions=30):
    report = {"tool": banner(), "date": datetime.date.today().isoformat(),
              "typer": _display(typer), "charge_engine": _display(charge_engine)}
    type_molecule = load_typer(typer) if {"B1", "B3"} & set(sets) else None
    if "B0" in sets:
        report["B0"] = b0.score_charges(b0.load(), load_charge_engine(charge_engine, toppar_dir(toppar)))
    if "QT" in sets:
        report["QT"] = training_check(toppar).as_dict()
    if "B1" in sets:
        cases, stats = b1.build(toppar_dir(toppar) / TOPOLOGY_FILE)
        report["B1"] = {"build": stats, "typing": score_typing(cases, type_molecule, n_confusions)}
    if "B3" in sets:
        cases, excluded = b3.build()
        report["B3"] = {"build": {"molecules": len(cases), "excluded": excluded},
                        "typing": score_typing(cases, type_molecule, n_confusions)}
    report["floor_violations"] = floors.check(report)
    return report


def _display(path):
    if path in (None, PRODUCTION):
        return PRODUCTION
    path = Path(path).resolve()
    return str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else str(path)


def to_markdown(report, confusions=True):
    out = [f"# Benchmark report", "", f"{report['tool']} · {report['date']}", "",
           f"- typer: `{report['typer']}`", f"- charge engine: `{report['charge_engine']}`", ""]
    if "B0" in report:
        r = report["B0"]
        out += ["## B0: charges vs CGenFF program 5.0 (P-v5 SI)", "",
                f"{r['exact']}/{r['molecules']} molecules exact (charges and penalties to 0.001)", "",
                "| molecule | atoms | analog params | max abs dq | max abs d penalty | exact |", "|---|---|---|---|---|---|"]
        out += [f"| {m['name']} | {m['atoms']} | {m['analog_parameters']} | {m['max_abs_dq']:.3f} | "
                f"{m['max_abs_dpenalty']:.3f} | {'yes' if m['exact'] else 'NO'} |" for m in r["per_molecule"]]
        out.append("")
    if "QT" in report:
        r = report["QT"]
        top = ", ".join(f"{x['type']} {x['q0']:+d} ({x['count']})" for x in r["formal_charges"][:12])
        out += ["## QT: increment charges vs training-residue charges (FINDINGS 2.2)", "",
                f"{r['residues']} residues with complete increments; median per-molecule RMSE {r['median_rmse']:.4f} e; "
                f"atoms within 0.01 e: {r['fraction_within_0.01']:.1%}, within 0.05 e: {r['fraction_within_0.05']:.1%}",
                "", f"Most frequent inferred formal charges: {top}", ""]
    for name in ("B1", "B3"):
        if name not in report:
            continue
        out += [f"## {name}: typing", "", f"build: `{json.dumps(report[name]['build'])}`", "",
                "| split | molecules fully correct | atoms correct | crashes | molecules/s |", "|---|---|---|---|---|"]
        for split, m in report[name]["typing"].items():
            out.append(f"| {split} | {m['molecule_accuracy']:.2%} ({m['molecules_correct']}/{m['molecules']}) | "
                       f"{m['atom_accuracy']:.2%} ({m['atoms_correct']}/{m['atoms']}) | {m['crashes']} | "
                       f"{m['molecules_per_second']:.0f} |")
        out.append("")
        for split, m in report[name]["typing"].items():
            if confusions and m.get("confusions"):
                out += [f"Top confusions, {name}-{split} (reference -> predicted):", "",
                        "| count | reference | predicted | example |", "|---|---|---|---|"]
                out += [f"| {c['count']} | {c['ref']} | {c['pred']} | `{c['examples'][0]['smiles']}` "
                        f"[{c['examples'][0]['name']}] |" for c in m["confusions"]]
                out.append("")
    v = report["floor_violations"]
    out += ["## Regression floors (SPEC 7.4)", "", "all met" if not v else "\n".join(f"- FAIL: {x}" for x in v), ""]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sets", default=",".join(DEFAULT_SETS), help="comma-separated benchmark sets")
    ap.add_argument("--typer", default=str(DEFAULT_TYPER), help="typer plug-in file defining type_molecule(mol)")
    ap.add_argument("--charge-engine", default=PRODUCTION,
                    help=f"'{PRODUCTION}' (default) or a module file such as prototype/charges.py")
    ap.add_argument("--toppar", help="CHARMM toppar directory (default: data/ from scripts/fetch_data.sh)")
    ap.add_argument("--out", default=str(REPO_ROOT / "build" / "report"), help="output directory")
    ap.add_argument("--confusions", type=int, default=30, help="confusion pairs listed per development split")
    ap.add_argument("--check-floors", action="store_true", help="exit 1 if any regression floor is violated")
    a = ap.parse_args(argv)
    report = run([s.strip() for s in a.sets.split(",")], a.typer, a.charge_engine, a.toppar, a.confusions)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=1))
    md = to_markdown(report)
    (out / "report.md").write_text(md)
    print(to_markdown(report, confusions=False))
    print(f"report: {out / 'report.md'}, {out / 'report.json'}")
    return 1 if a.check_floors and report["floor_violations"] else 0


if __name__ == "__main__":
    sys.exit(main())

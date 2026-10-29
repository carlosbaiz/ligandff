"""Benchmark runner: scores typing and charges on B0/B1/B2/B3 and writes JSON + markdown reports (SPEC 7.3).
Set QT is the charge formula check on the training residues (FINDINGS 2.2); set QF checks the production typing-rule
formal charges against the charges implied by the B1/B2 reference residues (SPEC 5.7); set QI checks the improper
centres from typing-rule `impr` actions against the IMPR lines of those residues (SPEC 5.13); set B4 is the
leave-one-out test of the analogy search on the program's own analogies recorded in the parameter file (SPEC 5.10);
set OUT (not run by default, needs OpenMM) builds an OpenMM system from the output of every B1/B2/B3 molecule.

usage: python -m open_cgenff.benchmarks [--sets B0,B1,B2,B3,B4,QT,QF,QI[,OUT]] [--typer FILE] [--charge-engine FILE]
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
from . import b0, b1, b2, b3, b4, floors, qf, qi
from . import out as output
from .common import load_typer
from .scoring import score_typing

PRODUCTION_TYPER = "open_cgenff.typing"
DEFAULT_TYPER = PRODUCTION_TYPER
DEFAULT_SETS = ("B0", "B1", "B2", "B3", "B4", "QT", "QF", "QI")
PRODUCTION = "open_cgenff.charges"
PILOT_TYPER = REPO_ROOT / "prototype" / "typer_pilot.py"            # frozen baseline (SPEC 7.4)


def resolve_typer(typer):
    """Production typer, or a plug-in file defining type_molecule(rdkit_mol) -> types."""
    if typer in (None, PRODUCTION_TYPER):
        from ..typing import type_rdkit
        return type_rdkit
    return load_typer(typer)


def _production_typed(rdmol):
    from ..molecule import InputMolecule
    from ..typing import type_molecule
    return type_molecule(InputMolecule.from_rdkit(rdmol))


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
    type_molecule = resolve_typer(typer) if {"B1", "B2", "B3"} & set(sets) else None
    built = {}
    if "B0" in sets:
        report["B0"] = b0.score_charges(b0.load(), load_charge_engine(charge_engine, toppar_dir(toppar)))
    if "QT" in sets:
        report["QT"] = training_check(toppar).as_dict()
    for name, module in (("B1", b1), ("B2", b2)):
        if {name, "QF", "QI", "OUT"} & set(sets):
            built[name] = module.build(toppar_dir(toppar) / TOPOLOGY_FILE)
        if name in sets:
            cases, stats = built[name]
            report[name] = {"build": stats, "typing": score_typing(cases, type_molecule, n_confusions)}
    if "QF" in sets:
        report["QF"] = qf.score_formal_charges({k: v[0] for k, v in built.items()}, _production_typed,
                                               IncrementTable.from_toppar(toppar_dir(toppar)), n_confusions)
    if "QI" in sets:
        report["QI"] = qi.score_impropers({k: v[0] for k, v in built.items()}, _production_typed, n_confusions)
    if "B3" in sets or "OUT" in sets:
        built["B3"] = b3.build()
    if "B3" in sets:
        cases, excluded = built["B3"]
        report["B3"] = {"build": {"molecules": len(cases), "excluded": excluded},
                        "typing": score_typing(cases, type_molecule, n_confusions)}
    if "B4" in sets:
        from ..params import default_engine
        engine = default_engine(toppar)
        report["B4"] = b4.score(engine, b4.build(engine.db), n_examples=n_confusions)
    if "OUT" in sets:
        report["OUT"] = output.score_output({k: v[0] for k, v in built.items()}, toppar)
    report["floor_violations"] = floors.check(report)
    return report


def _display(path):
    if path in (PRODUCTION, PRODUCTION_TYPER):
        return path
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
    for name in ("B1", "B2", "B3"):
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
    if "QF" in report:
        out += ["## QF: typing-rule formal charges vs reference charges (SPEC 5.7)", "",
                "Molecules with complete increments, no lone pairs and all types correct (production typer).", "",
                "| set | split | molecules | formal charges all correct | atoms wrong |", "|---|---|---|---|---|"]
        for name, r in report["QF"].items():
            for split, m in r["splits"].items():
                out.append(f"| {name} | {split} | {m['molecules']} | {m['correct']} | {m['atoms_wrong']} |")
        out.append("")
        rows = [(name, x) for name, r in report["QF"].items() for x in r["mismatches"]]
        if confusions and rows:
            out += ["Mismatches, development splits:", "", "| set | type | inferred | rules | count | example |",
                    "|---|---|---|---|---|---|"]
            out += [f"| {n} | {x['type']} | {x['inferred']:+d} | {x['rules']:+d} | {x['count']} | {x['example']} |"
                    for n, x in rows]
            out.append("")
    if "QI" in report:
        out += ["## QI: improper centres from typing rules vs residue IMPR lines (SPEC 5.13)", "",
                "Molecules with all types correct (production typer).", "",
                "| set | split | molecules | improper centres all correct | reference centres | missed | extra |",
                "|---|---|---|---|---|---|---|"]
        for name, r in report["QI"].items():
            for split, m in r["splits"].items():
                out.append(f"| {name} | {split} | {m['molecules']} | {m['correct']} | {m['reference_centres']} | "
                           f"{m['missed']} | {m['extra']} |")
        out.append("")
        rows = [(name, x) for name, r in report["QI"].items() for x in r["mismatches"]]
        if confusions and rows:
            out += ["Mismatches, development splits:", "", "| set | centre type | kind | count | example |",
                    "|---|---|---|---|---|"]
            out += [f"| {n} | {x['type']} | {x['kind']} | {x['count']} | {x['example']} |" for n, x in rows]
            out.append("")
    if "B4" in report:
        out += ["## B4: analogy search, leave-one-out on the program's recorded analogies (SPEC 5.10)", "",
                "| split | cases | same source | TPS within 0.5 | both |", "|---|---|---|---|---|"]
        for split, m in sorted(report["B4"]["splits"].items()):
            out.append(f"| {split} | {m['cases']} | {m['same_source_fraction']:.1%} | "
                       f"{m['tps_within_0.5_fraction']:.1%} | {m['both_fraction']:.1%} |")
        out.append("")
        if confusions and report["B4"]["examples"]:
            out += ["Failures, B4-train:", "", "| missing | program's source | found |", "|---|---|---|"]
            out += [f"| {x['case']} | {x['reference'][0]} ({x['reference'][1]}) | "
                    f"{x['found'][0] + ' (' + str(x['found'][1]) + ')' if x['found'] else '-'} |"
                    for x in report["B4"]["examples"]]
            out.append("")
    if "OUT" in report:
        r = report["OUT"]
        cols = ["molecules", "ok", "pipeline_errors", "openmm_built", "with_lone_pairs", "with_analogs",
                "missing_parameters"] + [f"param_penalty {b}" for _, b in output.PENALTY_BINS]
        out_md = ["## OUT: pipeline output builds an OpenMM system (SPEC 1, 5.14)", "",
                  f"pipeline (typing, parameters, charges, lone pairs): {r['molecules_per_second']:.1f} molecules/s per "
                  f"core; OpenMM build from `.str` + PSF + CGenFF files, {r['processes']} processes", "",
                  "| set/split | " + " | ".join(c.replace("_", " ") for c in cols) + " |", "|---" * (len(cols) + 1) + "|"]
        out_md += [f"| {k} | " + " | ".join(str(m.get(c, 0)) for c in cols) + " |" for k, m in r["splits"].items()]
        out_md.append("")
        if confusions and r["failures"]:
            out_md += ["Failures, development splits:", ""] + [f"- {x['set']} {x['name']}: {'; '.join(x['errors'])}"
                                                              for x in r["failures"]] + [""]
        out += out_md
    v = report["floor_violations"]
    out += ["## Regression floors (SPEC 7.4)", "", "all met" if not v else "\n".join(f"- FAIL: {x}" for x in v), ""]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sets", default=",".join(DEFAULT_SETS), help="comma-separated benchmark sets")
    ap.add_argument("--typer", default=DEFAULT_TYPER,
                    help=f"'{PRODUCTION_TYPER}' (default) or a plug-in file defining type_molecule(mol), "
                         f"e.g. prototype/typer_pilot.py")
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

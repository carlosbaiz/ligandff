# CLAUDE.md

Open-source, clean-room reimplementation of the CGenFF program (v5.0): atom typing, charges, bonded parameters by
analogy, CHARMM `.str` output. Read `docs/SPEC.md` before any work; it defines components, benchmarks, milestones
(M0–M7) and acceptance criteria. `docs/FINDINGS.md` has the evidence behind every number in the spec.

## Commands
```
bash scripts/fetch_data.sh                 # pinned CGenFF/CHARMM data into data/ (needed by tests)
python -m pytest -q                        # must stay green
python -m open_cgenff.benchmarks --check-floors   # B0/B1/B2/B3/QT/QF report (production typer) -> build/report/
python -m open_cgenff.benchmarks --sets B1,B2,B3 --typer prototype/typer_pilot.py   # pilot baseline
bash scripts/import_pi_si.sh <ci300363c_si_001.zip>  # P-I SI special cases -> data/pi_si/ (local only)
python prototype/build_rtf_benchmark.py data/toppar_c36_releases/charmm/toppar/toppar/top_all36_cgenff.rtf build/b1.pkl
python prototype/score.py build/b1.pkl prototype/typer_pilot.py --examples 30          # train-split failures only
python prototype/charmmgui_benchmark.py benchmarks/charmm_gui_cpa/raw prototype/typer_pilot.py
python prototype/charges.py fixtures --toppar data/toppar_c36_releases/charmm/toppar/toppar --fixtures tests/fixtures/si_v5_reference.json
```

## Hard rules
- Clean room: never download, run, decompile or inspect the official CGenFF program, its rule files or SilcsBio
  binaries (incl. `data/cgenff5/developer_program/*.exe`). Use only papers, SI and the public data in `data/`.
- Never commit `data/`, PDFs, or other third-party force-field files. Read force-field values from `data/` at runtime;
  never hand-edit parameter values or charge increments.
- Never call the tool "CGenFF"; it is "CGenFF-compatible". Outputs must say they are not from the official program.
- Benchmark hygiene: inspect per-molecule failures only on dev/train splits (B1-train, B2-train, B5-dev). For B1-test,
  B2-test and B5-test report aggregate metrics only. Never edit reference data to make a test pass.
- Every typing or penalty-rule change: rerun B0/B1/B2/B3, no unexplained regressions, log it in
  `docs/RULES_CHANGELOG.md` (rule, motivating molecules, metric delta).
- `prototype/` is the frozen pilot baseline: do not import it from `src/`; port code instead. Production must beat the
  regression floors in `docs/SPEC.md` §7.4.
- If a spec statement turns out wrong, fix the code, then update `docs/SPEC.md` (and FINDINGS) in the same change.

## Layout
`docs/` spec, findings, references, rules changelog · `prototype/` pilot code (reference) · `tests/` pytest +
fixtures · `benchmarks/` benchmark data and split manifests · `scripts/` data fetch · `src/open_cgenff/` production package (typing rules: `src/open_cgenff/typing/rules/cgenff_v5.rules`; penalty rules: `src/open_cgenff/params/rules/penalties_v5.rules`) · `.github/workflows/ci.yml` CI (tests + floors)

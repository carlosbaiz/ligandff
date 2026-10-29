# Prototype (pilot) code: frozen reference baseline

Written during the feasibility study (2026-10-29). Not production code: do not import from `src/`; port instead.
Ported in M0: `build_rtf_benchmark.py`, `score.py`, `charmmgui_benchmark.py` → `src/open_cgenff/benchmarks/`.
Bug fix (M0): `charmmgui_benchmark.py` finds topologies stored only as `cgen_<A>.itp` (git-ignored `<A>.itp`).
Paths below assume `TP=data/toppar_c36_releases/charmm/toppar/toppar` after `scripts/fetch_data.sh`.

| file | what it does | status / result |
|---|---|---|
| `charges.py` | charge increments + per-atom penalties (P-II); `fixtures` and `training` subcommands | exact on B0 (3 molecules, program 5.0); training check median 0.0098 e |
| `perceive.py` | rings ≤ 7, lowest-penalty resonance structure (via RDKit enumeration, not P-I's literal algorithm), ring classes/aromaticity (P-I) | used by the pilot typer |
| `typer_pilot.py` | ~110 typing rules in Python | B1 train 96.7 % / test 93.8 %; B3 34/34 |
| `build_rtf_benchmark.py` | builds B1 from `top_all36_cgenff.rtf` (+ split manifest) | 1,130 molecules |
| `score.py` | scores a typer module on B1; `--examples N` prints top train-split confusions | — |
| `query.py` | tabulates environments of given reference types on B1-train | rule-writing aid |
| `charmmgui_benchmark.py` | builds and scores B3 from CHARMM-GUI `.itp` files | 34 species |
| `penalty_annotations.py` | extracts single-substitution penalty-matrix elements from parameter-file annotations | 107 bonded + 280 nonbonded elements |

```
python prototype/charges.py fixtures --toppar $TP --fixtures tests/fixtures/si_v5_reference.json
python prototype/charges.py training --toppar $TP
python prototype/build_rtf_benchmark.py $TP/top_all36_cgenff.rtf build/b1.pkl
python prototype/score.py build/b1.pkl prototype/typer_pilot.py --examples 30
python prototype/charmmgui_benchmark.py benchmarks/charmm_gui_cpa/raw prototype/typer_pilot.py --detail
python prototype/penalty_annotations.py $TP/par_all36_cgenff_increments.prm
```

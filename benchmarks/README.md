# Benchmarks

Definitions, status and hygiene rules: `docs/SPEC.md` §7. Contents of this folder:

- `b1_rtf_zinc.split.json` – fixed train/test membership of B1 (built from data by `prototype/build_rtf_benchmark.py`;
  the molecules themselves come from the fetched RTF, so only names are stored here). B1-test is semi-sealed.
- `charmm_gui_cpa/` – B3: CHARMM-GUI GROMACS topologies supplied by Carlos Baiz (2026-10-29). See its `NOTES.md`.
- B2 (legacy RTF residues) and B5 (official CGenFF 5.0 outputs) will be added here (M2, M7). Planned B5 layout:
  `b5_official_v5/<name>.mol2` + `<name>.str`, plus `split.json` (dev 30 % / test 70 %, fixed seed).

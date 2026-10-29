# open-cgenff (working name)

An open-source, **CGenFF-compatible** atom typer and parameter/charge assigner for small molecules: a clean-room
reimplementation of the CGenFF program (v5.0) from the published algorithms and the public CGenFF 5.0 data.
Not affiliated with, or endorsed by, the CGenFF developers or SilcsBio. Status: specification, verified prototypes, package skeleton with benchmark runner (M0).

- `docs/SPEC.md` – technical specification, benchmarks, milestones
- `docs/FINDINGS.md` – what has been verified so far, with numbers
- `docs/REFERENCES.md` – papers and data sources

Quick start:
```
pip install -e ".[dev]"
bash scripts/fetch_data.sh      # fetches pinned CHARMM36/CGenFF 5.0 data into data/
python -m pytest -q
python -m open_cgenff.benchmarks --check-floors   # benchmark report in build/report/
```

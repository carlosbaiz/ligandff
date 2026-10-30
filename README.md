# LigandFF

An open-source, **CGenFF-compatible** atom typer and parameter/charge assigner for small molecules: a clean-room
reimplementation of the CGenFF program (v5.0) from the published algorithms and the public CGenFF 5.0 data.
Not affiliated with, or endorsed by, the CGenFF developers or SilcsBio. Status: specification, verified prototypes, package skeleton with benchmark runner (M0),
production charge engine (M1), perception and benchmark B2 (M2), rule-language typer with rule file v1 (M3),
parameter database and analogy engine (M4), penalty matrix fitted to the program's recorded analogies (M5), lone pairs,
`.str`/JSON/PSF output, CLI and OpenMM/GROMACS export (M6), final validation (M7, `docs/VALIDATION.md`).

- `docs/SPEC.md` – technical specification, benchmarks, milestones
- `docs/FINDINGS.md` – what has been verified so far, with numbers
- `docs/VALIDATION.md` – final validation (M7): section 1 targets, train/test, limits
- `docs/REFERENCES.md` – papers and data sources

Quick start:
```
pip install -e ".[dev]"
bash scripts/fetch_data.sh      # fetches pinned CHARMM36/CGenFF 5.0 data into data/
python -m pytest -q
python -m ligandff.benchmarks --check-floors   # benchmark report in build/report/
python -m ligandff.benchmarks.validation       # all sets + OpenMM build of every molecule; SPEC section 1 targets
```

Usage (input: mol2/SDF with all hydrogens and bond orders, or SMILES):
```
ligandff ligand.mol2 -o ligand.str --json ligand.json      # CHARMM stream file + report
ligandff "Clc1ccccc1C(=O)O" --name CBA -o cba.str --psf cba.psf
pip install -e ".[export]"; ligandff ligand.mol2 -o ligand.str --gromacs ligand   # .top/.gro via ParmEd
```
Penalties: < 10 fair analogy; 10–50 some validation; > 50 extensive validation. Outputs are CGenFF-compatible, not
from the official CGenFF program.

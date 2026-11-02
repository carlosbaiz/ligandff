# LigandFF

LigandFF is an open-source, **CGenFF-compatible** atom typer and parameter/charge assigner for small molecules. It is
a clean-room reimplementation of the CGenFF program (v5.0), written from the published algorithms and the public
CGenFF 5.0 data. For a molecule it assigns CGenFF atom types, partial charges with charge penalties, and bonded
parameters (by analogy where the force field has none, with parameter penalties). It writes a CHARMM stream file
(`.str`), a PSF, a JSON report, and OpenMM or GROMACS inputs.

LigandFF is **not** the official CGenFF program and is not affiliated with or endorsed by the CGenFF developers or
SilcsBio. Describe its output as "CGenFF-compatible parameters generated with LigandFF".

## Installation

You need Linux, macOS or Windows with WSL; `git`; and Python 3.10–3.13. Check with `python3 --version`. Run every
command from a terminal, in the order given.

**1. Download the code**
```bash
git clone https://github.com/carlosbaiz/ligandff.git
cd ligandff
```
Stay in the `ligandff` folder for the remaining steps.

**2. Create a Python environment and install LigandFF.** Use option A (plain Python) or option B (conda), not both.

Option A, plain Python:
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[export,dev]"
```

Option B, conda (Miniconda or Anaconda installed):
```bash
conda create -n ligandff python=3.12 -y
conda activate ligandff
python -m pip install -e ".[export,dev]"
```

This installs LigandFF with RDKit, NetworkX, OpenMM, ParmEd and pytest. Without `[export,dev]` you get only the core
(no OpenMM/GROMACS export, no tests).

**3. Download the CGenFF force-field files** (about 650 MB, a few minutes):
```bash
bash scripts/fetch_data.sh
```
LigandFF reads the CHARMM36/CGenFF 5.0 files from their public repositories at pinned versions. They are not
redistributed with LigandFF. The script puts them in `data/` inside the `ligandff` folder.

**4. Check the installation** (about one minute):
```bash
python -m pytest -q
```
The last line should report `passed` and no `failed`.

**5. Run a first molecule:**
```bash
ligandff "CC(=O)Nc1ccc(Cl)cc1" --name ACE -o ace.str --json ace.json --psf ace.psf
```
This writes `ace.str`, `ace.json` and `ace.psf` in the current folder.

**In a new terminal**, reactivate the environment before using LigandFF: `cd ligandff`, then
`source .venv/bin/activate` (option A) or `conda activate ligandff` (option B).

**From another folder:** the `ligandff` command finds the force-field files in the folder you installed it in. To use
another copy of the files, pass `--toppar /path/to/toppar` or set the environment variable `LIGANDFF_TOPPAR` to the
folder that contains `top_all36_cgenff.rtf`, `par_all36_cgenff.prm` and `par_all36_cgenff_increments.prm`.

**Update** to the latest version: `cd ligandff`, `git pull`, then step 4 again.

## Usage

Input: a mol2, SDF or MOL file with all hydrogens and bond orders, or a SMILES string (hydrogens are added).
```bash
ligandff ligand.mol2 -o ligand.str --json ligand.json     # CHARMM stream file + JSON report
ligandff ligand.mol2 -o ligand.str --psf ligand.psf       # + PSF (for CHARMM or OpenMM)
ligandff ligand.mol2 -o ligand.str --gromacs ligand       # + ligand.top / ligand.gro via ParmEd
ligandff --help                                           # all options
```
Use the `.str` file together with the CHARMM36/CGenFF files (`top_all36_cgenff.rtf`, `par_all36_cgenff.prm`), as
with any CGenFF stream file. GROMACS export does not support molecules with halogen lone pairs (Cl, Br, I on
aromatic rings); for those, convert the `.str` with the Lemkul lab's `cgenff_charmm2gmx`.

From Python, including an OpenMM system:
```python
from ligandff.pipeline import parameterize, read_molecule
from ligandff.export import openmm_system, potential_energy

result = parameterize(read_molecule("CC(=O)Nc1ccc(Cl)cc1", "ACE"), "ACE")
print(max(a.penalty for a in result.atoms))       # largest charge penalty
system, psf, positions = openmm_system(result)    # openmm.System from the CGenFF files + our .str
print(potential_energy(system, positions), "kJ/mol")
```

**Penalties** (as in CGenFF): below 10, the analogy is fair; 10–50, some validation is advised; above 50, extensive
validation is needed. Each parameter and charge in the `.str` and JSON outputs carries its penalty and, for
analogies, the source parameter it was copied from.

## Accuracy

On held-out molecules from the CGenFF training sets, LigandFF assigns every atom type correctly in 97–98 % of
molecules (99.6–99.7 % of atoms), reproduces the published CGenFF 5.0 charges and penalties exactly, and chooses the
same source parameter as the official program in 92 % of the analogies recorded in the parameter file. Every molecule
in the benchmarks builds a valid OpenMM system. Details and limitations: `docs/VALIDATION.md` and the technical report.

## Citation

If you use LigandFF, please cite the technical report (ChemRxiv, 2026; see `CITATION.cff` or GitHub's "Cite this
repository" button) together with the CGenFF papers:

- Vanommeslaeghe, K. et al. *J. Comput. Chem.* **2010**, 31, 671–690. doi:10.1002/jcc.21367
- Vanommeslaeghe, K.; MacKerell, A. D., Jr. *J. Chem. Inf. Model.* **2012**, 52, 3144–3154. doi:10.1021/ci300363c
- Vanommeslaeghe, K.; Raman, E. P.; MacKerell, A. D., Jr. *J. Chem. Inf. Model.* **2012**, 52, 3155–3168.
  doi:10.1021/ci3003649
- Croitoru, A. et al. *J. Chem. Theory Comput.* **2025**, 21, 3044–3065. doi:10.1021/acs.jctc.5c00046

## For developers

```bash
python -m ligandff.benchmarks --check-floors   # benchmark report (B0–B4) in build/report/
python -m ligandff.benchmarks.validation       # + OpenMM build of every benchmark molecule (~4 min)
```
- `docs/SPEC.md`: specification, benchmarks and acceptance criteria
- `docs/FINDINGS.md`: the evidence behind each number in the spec
- `docs/VALIDATION.md`: validation report (targets, train/test, limitations)
- `docs/RULES_CHANGELOG.md`: every typing and penalty-rule change with its benchmark effect
- `paper/`: source of the technical report and the scripts that regenerate its figures

## License

MIT (see `LICENSE`). The CHARMM/CGenFF force-field files are not part of LigandFF; `scripts/fetch_data.sh` downloads
them from their public repositories, and their own terms apply.

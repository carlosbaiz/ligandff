# LigandFF

LigandFF is an open-source, **CGenFF-compatible** atom typer and parameter/charge assigner for small molecules. It is
a clean-room reimplementation of the CGenFF program (v5.0), written from the published algorithms and the public
CGenFF 5.0 data. For a molecule it assigns CGenFF atom types, partial charges with charge penalties, and bonded
parameters (by analogy where the force field has none, with parameter penalties). It writes a CHARMM stream file
(`.str`), a PSF, a JSON report, and OpenMM or GROMACS inputs.

LigandFF is **not** the official CGenFF program and is not affiliated with or endorsed by the CGenFF developers or
SilcsBio. Describe its output as "CGenFF-compatible parameters generated with LigandFF".

## Installation

You need Linux, macOS or Windows with WSL (Windows Subsystem for Linux), and `git`. Run every command below in a
terminal, in the order given. There are two ways to install: option A uses conda and is recommended; option B uses
plain Python 3.10–3.13. Both give the same program.

### Option A: conda (recommended)

You need Miniconda, Miniforge or Anaconda (check with `conda --version`).

**1. Download the code**
```bash
git clone https://github.com/carlosbaiz/ligandff.git
cd ligandff
```
Stay in the `ligandff` folder for the remaining steps.

**2. Create an environment and install LigandFF**
```bash
conda create -n ligandff --override-channels -c conda-forge python=3.12 rdkit networkx openmm pytest -y
conda activate ligandff
python -m pip install -e ".[export,dev]"
```
`--override-channels` makes conda use only the free conda-forge channel. Without it, recent Miniconda and Anaconda
installations also query Anaconda's default channels and stop with `CondaToSNonInteractiveError` until their Terms of
Service are accepted.

**3. Download the CGenFF force-field files** (about 650 MB, under a minute on a fast connection)
```bash
bash scripts/fetch_data.sh
```

**4. Check the installation** (about one minute)
```bash
python -m pytest -q
```
The last line should read like `124 passed, 3 skipped` and must not mention `failed`. Two skipped tests need local
data that is not distributed; the third runs only where GROMACS itself (`gmx`) is installed.

**5. Run a first molecule**
```bash
ligandff "CC(=O)Nc1ccc(Cl)cc1" --name LIG -o lig.str --json lig.json --psf lig.psf
```
This prints `LIG: param penalty 0.0, charge penalty 0.0, 0 parameters by analogy` and writes, in the current
folder, `lig.str` (CHARMM), `lig.json` (report), `lig.psf` (CHARMM/OpenMM) and the GROMACS files `lig.itp`, `lig.top`
and `lig.gro`.

In every new terminal, run `conda activate ligandff` before using LigandFF.

### Option B: plain Python (no conda)

You need Python 3.10–3.13 (check with `python3 --version`).
```bash
git clone https://github.com/carlosbaiz/ligandff.git
cd ligandff
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[export,dev]"
bash scripts/fetch_data.sh
python -m pytest -q
ligandff "CC(=O)Nc1ccc(Cl)cc1" --name LIG -o lig.str --json lig.json --psf lig.psf
```
The test summary and the output are the same as in option A. In every new terminal, run `cd ligandff` and
`source .venv/bin/activate` before using LigandFF.

### Notes

- **Force-field files.** LigandFF reads the CHARMM36/CGenFF 5.0 files (`top_all36_cgenff.rtf`,
  `par_all36_cgenff.prm`, `par_all36_cgenff_increments.prm`). They are not redistributed with LigandFF;
  `scripts/fetch_data.sh` downloads them from their public repositories at pinned versions into `data/`. The
  `ligandff` command finds them there from any folder. To use another copy, pass `--toppar /path/to/folder` or set
  the environment variable `LIGANDFF_TOPPAR` to the folder that contains the three files.
- **Updating.** `cd ligandff`, `git pull`, then `python -m pytest -q` again.

## Usage

Input: a mol2, SDF or MOL file with all hydrogens and bond orders, or a SMILES string (hydrogens are added).
```bash
ligandff ligand.mol2 -o ligand.str                        # ligand.str + GROMACS ligand.itp, ligand.top, ligand.gro
ligandff ligand.mol2 -o ligand.str --json ligand.json     # + JSON report
ligandff ligand.mol2 -o ligand.str --psf ligand.psf       # + PSF (for CHARMM or OpenMM)
ligandff ligand.mol2 -o ligand.str --gromacs gmx/lig      # GROMACS files as gmx/lig.itp, .top, .gro instead
ligandff ligand.mol2 -o ligand.str --no-gromacs           # no GROMACS files
ligandff --help                                           # all options
```
**CHARMM / OpenMM:** use the `.str` file together with the CHARMM36/CGenFF files (`top_all36_cgenff.rtf`,
`par_all36_cgenff.prm`), as with any CGenFF stream file.

**GROMACS:** three files are written.
- `ligand.itp`: the molecule (`[ moleculetype ]`) with every bonded and 1-4 parameter written out explicitly.
  Halogen lone pairs (Cl, Br, I on aromatic rings) are virtual sites.
- `ligand.top`: a standalone topology (CGenFF atom types and NBFIX pairs, then `#include "ligand.itp"`) for the
  ligand alone, e.g. in vacuum or as the start of a solvated box.
- `ligand.gro`: starting coordinates, from the input file or generated by RDKit for SMILES.

For a protein–ligand system with a CHARMM36 force field for GROMACS, include `ligand.itp` after the force field
instead of using `ligand.top`. The force field must define the CGenFF atom types used in the `.itp` (and `LPH`
for molecules with halogen lone pairs). Check the `[ atomtypes ]` of `ligand.top` against it.

From Python, including an OpenMM system:
```python
from ligandff.pipeline import parameterize, read_molecule
from ligandff.export import openmm_system, potential_energy

result = parameterize(read_molecule("CC(=O)Nc1ccc(Cl)cc1", "LIG"), "LIG")
print(result.param_penalty, result.charge_penalty, len(result.analogs))
for atom in result.atoms:
    print(atom.name, atom.type, atom.charge, atom.charge_penalty)
system, psf, positions = openmm_system(result)    # openmm.System from the CGenFF files + our .str
print(potential_energy(system, positions), "kJ/mol")
```

**Penalties** (as in CGenFF): below 10, the analogy is fair; 10–50, some validation is advised; above 50, extensive
validation is needed. The `.str` file and the JSON report list each charge and parameter with its penalty and, for
parameters taken by analogy, the source parameter they were copied from.

Penalties marked **estimated** (`estimated penalty=` in the `.str`, `"estimated": true` in the JSON, "(estimated)" in
the command-line summary) are not on the official scale. The analogy needs a substitution never seen in CGenFF's
recorded analogies, and LigandFF can only guess its cost. Do not apply the thresholds above to them: validate or
optimize those terms whatever the number says.

## Accuracy

On held-out molecules from the CGenFF training sets, LigandFF assigns every atom type correctly in 97–98 % of
molecules (99.6–99.7 % of atoms), reproduces the published CGenFF 5.0 charges and penalties exactly, and, replaying the analogies recorded in the
parameter file under the conditions they were made, chooses the same source parameter as the official program in
92 % of cases (36 % for the few annotations made outside that batch; agreement on new molecules is not measured). Every molecule
in the benchmarks builds a valid OpenMM system, and its GROMACS files reproduce the same energies. Details and limitations: `docs/VALIDATION.md` and the technical report.

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

## License

MIT (see `LICENSE`). The CHARMM/CGenFF force-field files are not part of LigandFF; `scripts/fetch_data.sh` downloads
them from their public repositories, and their own terms apply.

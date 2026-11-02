"""OpenMM and GROMACS export (SPEC 5.14), reusing established readers/converters rather than new force-field code:
OpenMM's CHARMM PSF/parameter readers build the system from the CGenFF files plus our .str; ParmEd's (optional) CHARMM
readers do the same for GROMACS .top/.gro.
"""
import tempfile
from pathlib import Path

from .io.psf import write_psf
from .io.stream import write_str
from .params import default_engine
from .paths import PARAMETER_FILE, TOPOLOGY_FILE, toppar_dir


def coordinates(result, toppar=None, seed=42):
    """Input coordinates, or an RDKit embedding (Angstrom; a starting geometry only) for the real atoms: of the input
    molecule with its bond orders and formal charges when the result has one, else of the connectivity alone."""
    if result.coords is not None:
        return [tuple(c) for c in result.coords]
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mols = []
    if result.source is not None:
        mols.append(_rdkit_input(result.source))
    elements = default_engine(toppar).db.elements
    mols.append(_rdkit_graph([elements[a.type][1] for a in result.atoms], [(i, j, 1) for i, j in result.bonds]))
    for mol in mols:
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
            if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True) == 0:
                return [tuple(p) for p in mol.GetConformer().GetPositions()]
        except (ValueError, RuntimeError):            # RDKit valence errors (e.g. P- with four single bonds)
            continue
    raise ValueError(f"{result.name}: could not generate coordinates")


def _rdkit_graph(elements, bonds, charges=None):
    from rdkit import Chem
    from .molecule import AROMATIC
    kinds = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE, AROMATIC: Chem.BondType.AROMATIC}
    rw = Chem.RWMol()
    for k, el in enumerate(elements):
        atom = Chem.Atom(el)
        atom.SetNoImplicit(True)                          # every H is explicit
        atom.SetFormalCharge(charges[k] if charges else 0)
        rw.AddAtom(atom)
    for i, j, order in bonds:
        rw.AddBond(i, j, kinds[order])
        if order == AROMATIC:
            for k in (i, j):
                rw.GetAtomWithIdx(k).SetIsAromatic(True)
            rw.GetBondBetweenAtoms(i, j).SetIsAromatic(True)
    return rw.GetMol()


def _rdkit_input(mol):
    return _rdkit_graph(mol.elements, mol.bonds, mol.formal_charges)


def openmm_system(result, toppar=None, workdir=None, **create_args):
    """(openmm.System, CharmmPsfFile, positions in nm incl. lone pairs) for a pipeline Result.

    Writes <name>.psf and <name>.str to workdir (a temporary directory by default) and reads them with
    openmm.app.CharmmPsfFile / CharmmParameterSet(top_all36_cgenff.rtf, par_all36_cgenff.prm, <name>.str).
    create_args go to CharmmPsfFile.createSystem (default: no cutoff, no constraints).
    """
    import openmm
    from openmm import app, unit
    masses = default_engine(toppar).db.file.masses
    xyz = coordinates(result, toppar) + [(0.0, 0.0, 0.0)] * len(result.lone_pairs)
    toppar = toppar_dir(toppar)
    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory()
        workdir = tmp.name
    workdir = Path(workdir)
    psf_path, str_path = workdir / f"{result.name}.psf", workdir / f"{result.name}.str"
    write_psf(result, masses, psf_path)
    write_str(result, str_path)
    params = app.CharmmParameterSet(*_parameter_files(toppar, str_path))
    psf = app.CharmmPsfFile(str(psf_path))
    args = {"nonbondedMethod": app.NoCutoff, "constraints": None}
    args.update(create_args)
    system = psf.createSystem(params, **args)
    positions = unit.Quantity([openmm.Vec3(*p) * 0.1 for p in xyz], unit.nanometer)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    context.computeVirtualSites()
    positions = context.getState(getPositions=True).getPositions()
    if tmp is not None:
        tmp.cleanup()
    return system, psf, positions


def _parameter_files(toppar, str_path):
    return str(toppar / TOPOLOGY_FILE), str(toppar / PARAMETER_FILE), str(str_path)


def potential_energy(system, positions):
    """Potential energy (kJ/mol) on the Reference platform."""
    import openmm
    from openmm import unit
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def write_gromacs(result, prefix, toppar=None):
    """<prefix>.top and <prefix>.gro via ParmEd (optional dependency) from our PSF + the CGenFF files + our
    .str (ParmEd's CHARMM readers, so the CGenFF atom type names are kept).

    Not for molecules with halogen lone pairs: ParmEd does not carry CHARMM lone pairs into GROMACS virtual sites, so
    the lone pair would become a free particle. Use the .str with the Lemkul lab's cgenff_charmm2gmx converter for those."""
    if result.lone_pairs:
        raise NotImplementedError(f"{result.name}: GROMACS export of halogen lone pairs is not supported through "
                                  "ParmEd; convert the .str with cgenff_charmm2gmx instead")
    try:
        import parmed
    except ImportError as exc:
        raise RuntimeError("GROMACS export needs ParmEd: conda install -c conda-forge parmed (or pip install parmed, "
                           "which needs a C++ compiler)") from exc
    from parmed.charmm import CharmmParameterSet, CharmmPsfFile
    masses = default_engine(toppar).db.file.masses
    xyz = coordinates(result, toppar)
    with tempfile.TemporaryDirectory() as workdir:
        psf_path, str_path = Path(workdir) / f"{result.name}.psf", Path(workdir) / f"{result.name}.str"
        write_psf(result, masses, psf_path)
        write_str(result, str_path)
        structure = CharmmPsfFile(str(psf_path))
        structure.load_parameters(CharmmParameterSet(*_parameter_files(toppar_dir(toppar), str_path)))
    structure.coordinates = xyz                           # Angstrom; types and CHARMM terms kept from the PSF/files
    structure.save(f"{prefix}.top", format="gromacs", overwrite=True, combine="all")
    structure.save(f"{prefix}.gro", overwrite=True)
    return f"{prefix}.top", f"{prefix}.gro"

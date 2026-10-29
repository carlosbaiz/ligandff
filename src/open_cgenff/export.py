"""OpenMM and GROMACS export (SPEC 5.14), reusing established readers/converters rather than new force-field code:
OpenMM's CHARMM PSF/parameter readers build the system from the CGenFF files plus our .str; ParmEd (optional) converts
that system to GROMACS .top/.gro.
"""
import tempfile
from pathlib import Path

from .io.psf import write_psf
from .io.stream import write_str
from .params import default_engine
from .paths import PARAMETER_FILE, TOPOLOGY_FILE, toppar_dir


def coordinates(result, toppar=None, seed=42):
    """Input coordinates, or an RDKit embedding of the connectivity (Angstrom; a starting geometry only: bond orders
    are not used), for the real atoms."""
    if result.coords is not None:
        return [tuple(c) for c in result.coords]
    from rdkit import Chem
    from rdkit.Chem import AllChem
    elements = default_engine(toppar).db.elements
    rw = Chem.RWMol()
    for a in result.atoms:
        atom = Chem.Atom(elements[a.type][1])
        atom.SetNoImplicit(True)                          # every H is explicit
        rw.AddAtom(atom)
    for i, j in result.bonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True) != 0:
        raise ValueError(f"{result.name}: could not generate coordinates")
    return [tuple(p) for p in mol.GetConformer().GetPositions()]


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
    params = app.CharmmParameterSet(str(toppar / TOPOLOGY_FILE), str(toppar / PARAMETER_FILE), str(str_path))
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


def potential_energy(system, positions):
    """Potential energy (kJ/mol) on the Reference platform."""
    import openmm
    from openmm import unit
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def write_gromacs(result, prefix, toppar=None):
    """<prefix>.top and <prefix>.gro via ParmEd (optional dependency) from the OpenMM system."""
    try:
        import parmed
    except ImportError as exc:
        raise RuntimeError("GROMACS export needs ParmEd: pip install parmed") from exc
    system, psf, positions = openmm_system(result, toppar)
    structure = parmed.openmm.load_topology(psf.topology, system, xyz=positions)
    structure.save(f"{prefix}.top", overwrite=True)
    structure.save(f"{prefix}.gro", overwrite=True)
    return f"{prefix}.top", f"{prefix}.gro"

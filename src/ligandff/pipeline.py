"""End-to-end pipeline (SPEC section 4): read -> perceive -> type -> bonded parameters (exact or by analogy) ->
charges (with analog increments and penalties) -> halogen lone pairs.

    result = parameterize(read_molecule("ligand.mol2"))
    result.atoms, result.param_penalty, result.charge_penalty, result.analogs
"""
from dataclasses import dataclass, field
from pathlib import Path

from .charges import IncrementTable, charge_molecule
from .lonepairs import apply_lone_pair_charges, lone_pairs
from .molecule import InputMolecule
from .params import assign_parameters, charge_analogs, default_engine
from .typing import type_molecule

DEFAULT_RESNAME = "LIG"
_TABLES = {}


@dataclass
class Atom:
    name: str
    type: str
    charge: float
    charge_penalty: float
    formal_charge: int = 0
    charge_penalty_estimated: bool = False    # includes an analog penalty from unobserved substitutions


@dataclass
class Result:
    name: str                         # residue name
    atoms: list                       # [Atom], input order (lone pairs not included)
    bonds: list                       # [(i, j)]
    assignment: object                # params.ParameterAssignment (bonds, angles, dihedrals, impropers)
    lone_pairs: list                  # [lonepairs.LonePair]
    coords: list = None               # input coordinates (Angstrom), if any
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)       # typing errors: the output should not be used as is
    source: object = None             # the InputMolecule (bond orders for coordinate generation), if known

    @property
    def net_charge(self):
        return sum(a.formal_charge for a in self.atoms)

    @property
    def param_penalty(self):
        return self.assignment.param_penalty

    @property
    def charge_penalty(self):
        return max((a.charge_penalty for a in self.atoms), default=0.0)

    @property
    def param_penalty_estimated(self):
        """The RESI param penalty comes from an analog priced with unobserved substitutions (not the official scale)."""
        top = self.param_penalty
        return any(p.analog is not None and p.analog.estimated and p.penalty == top for p in self.assignment.parameters)

    @property
    def charge_penalty_estimated(self):
        top = self.charge_penalty
        return any(a.charge_penalty_estimated and a.charge_penalty == top for a in self.atoms)

    @property
    def analogs(self):
        """{(kind, types): AnalogResult}: parameters absent from CGenFF, in first-use order."""
        return self.assignment.analogs()

    @property
    def impropers(self):
        return [p for p in self.assignment.parameters if p.kind == "improper"]


def increment_table(toppar=None):
    if toppar not in _TABLES:
        _TABLES[toppar] = IncrementTable.from_toppar(toppar)
    return _TABLES[toppar]


def parameterize(mol, name=None, toppar=None, engine=None, table=None):
    """Full assignment for an InputMolecule (explicit hydrogens, bond orders or aromatic flags)."""
    words = (name or mol.name or "").split()
    result = parameterize_typed(type_molecule(mol), words[0][:6] if words else DEFAULT_RESNAME, toppar, engine, table,
                                mol.coords)
    result.source = mol
    return result


def parameterize_typed(typed, name=DEFAULT_RESNAME, toppar=None, engine=None, table=None, coords=None):
    """Parameters, charges and lone pairs for an already typed molecule (TypedMolecule)."""
    engine = engine or default_engine(toppar)
    table = table or increment_table(toppar)
    assignment = assign_parameters(typed, engine)
    charges = charge_molecule(typed, table, charge_analogs(typed, engine, table))
    lps = lone_pairs(typed.names, typed.types, typed.bonds)
    q = apply_lone_pair_charges(charges.charges, lps)
    atoms = [Atom(n, t, round(c, 6), round(p, 6), f, e) for n, t, c, p, f, e in
             zip(typed.names, typed.types, q, charges.penalties, typed.formal_charges, charges.estimated)]
    warnings = list(typed.warnings)
    warnings += [f"no charge increments for {' '.join(k)}" for k in charges.missing]
    warnings += [f"no parameter or analog for {kind} {' '.join(t)}" for kind, t in assignment.missing]
    return Result(name, atoms, list(typed.bonds), assignment, lps, coords, warnings, list(typed.errors))


def read_molecule(source, name=None):
    """InputMolecule from a .mol2 / .sdf / .mol file or a SMILES string (hydrogens added for SMILES; files must
    contain all hydrogens, as for the official program)."""
    from rdkit import Chem
    path = Path(source)
    if path.suffix.lower() == ".mol2" and path.exists():
        from .io.mol2 import read_mol2
        mols = read_mol2(path)
        if not mols:
            raise ValueError(f"{source}: no molecule")
        mol = mols[0]
        if name:
            mol.name = name
        return mol
    if path.suffix.lower() in (".sdf", ".mol") and path.exists():
        rd = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)), None)
        if rd is None:
            raise ValueError(f"{source}: could not read molecule")
        return InputMolecule.from_rdkit(rd, name or (rd.GetProp("_Name") if rd.HasProp("_Name") else path.stem))
    rd = Chem.MolFromSmiles(str(source))
    if rd is None:
        raise ValueError(f"{source}: not a readable file or SMILES")
    return InputMolecule.from_rdkit(Chem.AddHs(rd), name or "")

"""Shared benchmark definitions: cases, type equivalence, typer plug-ins."""
import importlib.util
from dataclasses import dataclass, field
from pathlib import Path

ELEMENT_PREFIXES = [("HG", "H"), ("CG", "C"), ("NG", "N"), ("OG", "O"), ("SG", "S"), ("PG", "P"), ("FG", "F"),
                    ("CLG", "Cl"), ("BRG", "Br"), ("IG", "I"), ("BG", "B"), ("SEG", "Se"), ("ALG", "Al")]

# 'altnum' types: the 1/2 digit is an arbitrary, parameter-symmetric choice (SPEC 5.5, 7.3)
ALTNUM = {"CG2DC1": "CG2DC?", "CG2DC2": "CG2DC?", "CG2D1O": "CG2D?O", "CG2D2O": "CG2D?O",
          "CG25C1": "CG25C?", "CG25C2": "CG25C?", "CG251O": "CG25?O", "CG252O": "CG25?O"}


def canonical_type(t):
    return ALTNUM.get(t, t)


def element_of_type(t):
    """Element of a CGenFF atom type from its prefix, or None (lone pairs, ions, non-CGenFF types)."""
    for prefix, el in ELEMENT_PREFIXES:
        if t.startswith(prefix):
            return el
    return None


@dataclass
class TypingCase:
    """One reference molecule: RDKit molecule with explicit H and reference types in RDKit atom order."""
    name: str
    mol: object
    types: list
    split: str = "dev"              # 'train' / 'test' (sealed or semi-sealed) / 'dev'
    smiles: str = None
    meta: dict = field(default_factory=dict)


def load_typer(path):
    """Load a typer plug-in from a file path; it must define type_molecule(rdkit_mol) -> list of types."""
    path = Path(path).resolve()
    spec = importlib.util.spec_from_file_location(f"typer_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.type_molecule

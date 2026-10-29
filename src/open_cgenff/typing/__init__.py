"""Atom typing (SPEC 5.5-5.7): the P-I rule language and the rule file.

    from open_cgenff.typing import type_molecule
    tm = type_molecule(InputMolecule.from_rdkit(mol))    # TypedMolecule: types, typing-rule formal charges, impropers
"""
from functools import lru_cache
from pathlib import Path

from ..molecule import TypedMolecule
from ..perception import perceive
from .engine import TypingError, TypingResult, type_perception
from .rules import RuleSet, RuleSyntaxError, parse

RULES_FILE = Path(__file__).parent / "rules" / "cgenff_v5.rules"


@lru_cache(maxsize=None)
def default_rules():
    return RuleSet.from_file(RULES_FILE)


def type_molecule(mol, rules=None):
    """Perceive and type an InputMolecule; returns a TypedMolecule (perception and rule warnings included)."""
    p = perceive(mol)
    res = type_perception(rules or default_rules(), p)
    bonds = [(i, j) for i, j, _ in mol.bonds]
    return TypedMolecule(list(mol.atom_names), res.types, bonds, res.charges, res.impropers,
                         p.warnings + res.warnings, res.errors)


def type_rdkit(rdmol):
    """Benchmark typer interface: RDKit molecule with explicit H -> list of types."""
    from ..molecule import InputMolecule
    return type_molecule(InputMolecule.from_rdkit(rdmol)).types


__all__ = ["RuleSet", "RuleSyntaxError", "TypingError", "TypingResult", "default_rules", "parse", "type_molecule",
           "type_perception", "type_rdkit", "RULES_FILE"]

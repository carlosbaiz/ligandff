"""Charge engine (SPEC 5.11): bond-charge increments and per-atom charge penalties.

    table = IncrementTable.from_toppar()
    result = charge_molecule(typed_molecule, table, analogs={missing_types: Analog(source_types, penalty)})
    result.charges, result.penalties, result.charge_penalty, result.missing
"""
from .increments import (Analog, ChargeResult, IncrementTable, assign_charges, charge_molecule,
                         parameter_instances, MAX_DIHEDRAL_PENALTY)
from .penalties import atom_penalty
from .training import TrainingCheck, training_check

__all__ = ["Analog", "ChargeResult", "IncrementTable", "assign_charges", "charge_molecule", "parameter_instances",
           "atom_penalty", "training_check", "TrainingCheck", "MAX_DIHEDRAL_PENALTY"]

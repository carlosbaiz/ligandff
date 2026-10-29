"""Parameter database, penalty model and analogy engine (SPEC 5.8-5.10).

    db = ParameterDatabase.from_toppar()
    engine = AnalogyEngine(db, fitted_model(db))
    assignment = assign_parameters(typed_molecule, engine)       # .parameters, .param_penalty, .analogs()
    analogs = charge_analogs(typed_molecule, engine, increment_table)   # for charges.charge_molecule
"""
from .analogy import (AnalogResult, AnalogyEngine, AssignedParameter, ParameterAssignment, assign_parameters,
                      charge_analogs, molecule_instances)
from .database import ParameterDatabase
from .penalties import (PENALTY_RULES_FILE, BondGroups, PenaltyModel, PenaltyRuleError, PenaltyRules,
                        SubstitutionMatrix, annotation_split, fitted_model, parse_penalty_rules,
                        provisional_model)

__all__ = ["AnalogResult", "AnalogyEngine", "AssignedParameter", "ParameterAssignment", "assign_parameters",
           "charge_analogs", "molecule_instances", "ParameterDatabase", "PENALTY_RULES_FILE", "BondGroups",
           "PenaltyModel", "PenaltyRuleError", "PenaltyRules", "SubstitutionMatrix", "parse_penalty_rules",
           "provisional_model", "fitted_model", "annotation_split", "default_engine"]

_DEFAULT = {}


def default_engine(toppar=None):
    """AnalogyEngine on the fetched parameter file with the fitted penalty model (M5; cached per directory)."""
    if toppar not in _DEFAULT:
        db = ParameterDatabase.from_toppar(toppar)
        _DEFAULT[toppar] = AnalogyEngine(db, fitted_model(db))
    return _DEFAULT[toppar]

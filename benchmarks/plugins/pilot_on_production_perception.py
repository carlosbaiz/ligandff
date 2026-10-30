"""Typer plug-in: the frozen pilot typer running on production perception (ligandff.perception).
SPEC 5.3 acceptance (pilot regression on B1 not worse). Glue outside src/, which never imports prototype/.
usage: python -m ligandff.benchmarks --typer benchmarks/plugins/pilot_on_production_perception.py"""
import importlib.util, os

from ligandff.molecule import InputMolecule
from ligandff.perception import perceive

PROTO = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'prototype')


def _load(name):
    spec = importlib.util.spec_from_file_location(f'adapter_{name}', os.path.join(PROTO, f'{name}.py'))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def production_perceive_factory():
    proto = _load('perceive')

    def perceive_production(rdmol):
        m = proto.perceive(rdmol)                       # keeps pilot-only helpers: nH, deg, nbrs, posany/negany
        p = perceive(InputMolecule.from_rdkit(rdmol))
        m.rings = p.rings
        m.bo = [{j: p.order(i, j) for j in p.neighbors[i]} for i in range(p.n_atoms)]
        m.charge = list(p.charges)
        m.rcls = dict(p.ring_class)
        m.val = [p.valence(i) for i in range(p.n_atoms)]
        m.atom_rings = [list(r) for r in p.atom_rings]
        m.inring = p.in_ring
        m.fused = [any(len(set(r1) & set(r2)) >= 2 and sum(j in set(r1) | set(r2) for j in m.nbrs[i]) >= 3
                       for a, r1 in enumerate(m.atom_rings[i]) for r2 in m.atom_rings[i][a + 1:]) for i in range(m.n)]
        m.perception = p
        return m
    return perceive_production


def pilot_typer_on_production_perception():
    typer = _load('typer_pilot')
    typer.perceive = production_perceive_factory()
    return typer.type_molecule


type_molecule = pilot_typer_on_production_perception()

"""Charge engine must reproduce CGenFF program 5.0 output exactly (charges and per-atom charge penalties)."""
import json, os
import charges


def test_si_reference_outputs_exact(toppar, root):
    table = charges.load_increments(os.path.join(toppar, 'par_all36_cgenff_increments.prm'))
    for mol in json.load(open(os.path.join(root, 'tests', 'fixtures', 'si_v5_reference.json'))):
        names = [a['name'] for a in mol['atoms']]; idx = {n: i for i, n in enumerate(names)}
        bonds = [(idx[a], idx[b]) for a, b in mol['bonds']]
        analogs = {tuple(x['missing']): (x['source'], x['penalty']) for x in mol['analog_parameters']}
        q, pen, missing = charges.assign([a['type'] for a in mol['atoms']], bonds,
                                         [mol['formal_charges'].get(n, 0) for n in names], table, analogs)
        assert not missing, (mol['name'], missing)
        for a, qi, pi in zip(mol['atoms'], q, pen):
            assert round(qi, 3) == a['charge'], (mol['name'], a['name'], qi, a['charge'])
            assert round(pi, 3) == a['charge_penalty'], (mol['name'], a['name'], pi, a['charge_penalty'])
        assert round(max(pen), 3) == mol['charge_penalty']

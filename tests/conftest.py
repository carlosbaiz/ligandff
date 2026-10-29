import os, sys, pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'prototype'))
TOPPAR = os.path.join(ROOT, 'data', 'toppar_c36_releases', 'charmm', 'toppar', 'toppar')


@pytest.fixture(scope='session')
def toppar():
    if not os.path.exists(os.path.join(TOPPAR, 'par_all36_cgenff_increments.prm')):
        pytest.skip('CGenFF data not fetched: run scripts/fetch_data.sh')
    return TOPPAR


@pytest.fixture(scope='session')
def root():
    return ROOT

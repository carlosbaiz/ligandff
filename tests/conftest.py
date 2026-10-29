import os, pickle, subprocess, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'prototype'))
sys.path.insert(0, os.path.join(ROOT, 'src'))           # production package without installation
TOPPAR = os.path.join(ROOT, 'data', 'toppar_c36_releases', 'charmm', 'toppar', 'toppar')


@pytest.fixture(scope='session')
def toppar():
    if not os.path.exists(os.path.join(TOPPAR, 'par_all36_cgenff_increments.prm')):
        if os.environ.get('OPEN_CGENFF_REQUIRE_DATA'):    # set in CI: missing data must fail, not skip
            pytest.fail('CGenFF data not fetched: run scripts/fetch_data.sh')
        pytest.skip('CGenFF data not fetched: run scripts/fetch_data.sh')
    return TOPPAR


@pytest.fixture(scope='session')
def root():
    return ROOT


@pytest.fixture(scope='session')
def b1(toppar, root, tmp_path_factory):
    """B1 as built by the frozen prototype builder."""
    out = tmp_path_factory.mktemp('b1') / 'b1.pkl'
    subprocess.run([sys.executable, os.path.join(root, 'prototype', 'build_rtf_benchmark.py'),
                    os.path.join(toppar, 'top_all36_cgenff.rtf'), str(out)], check=True, capture_output=True)
    return pickle.load(open(out, 'rb'))

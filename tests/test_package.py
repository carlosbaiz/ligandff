"""Package skeleton (M0): importable, versioned, and every output banner carries the SPEC 2.2 disclaimer."""
import open_cgenff
from open_cgenff.paths import toppar_dir


def test_banner_disclaims_official_program():
    b = open_cgenff.banner()
    assert open_cgenff.__version__ in b and 'CGenFF-compatible' in b and 'not the official CGenFF program' in b
    assert 'CGenFF version 5.0' in b


def test_toppar_dir_finds_fetched_data(toppar):
    assert str(toppar_dir()) == toppar

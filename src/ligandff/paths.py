"""Locations of the fetched third-party data (scripts/fetch_data.sh) and of repository-relative files."""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TOPPAR_ENV = "LIGANDFF_TOPPAR"
INCREMENTS_FILE = "par_all36_cgenff_increments.prm"
TOPOLOGY_FILE = "top_all36_cgenff.rtf"
PARAMETER_FILE = "par_all36_cgenff.prm"


class DataNotFound(FileNotFoundError):
    pass


def toppar_dir(path=None):
    """CHARMM toppar directory: explicit argument, else $LIGANDFF_TOPPAR, else <repo>/data/... from fetch_data.sh."""
    cand = Path(path or os.environ.get(TOPPAR_ENV) or
                REPO_ROOT / "data" / "toppar_c36_releases" / "charmm" / "toppar" / "toppar")
    if not (cand / INCREMENTS_FILE).exists():
        raise DataNotFound(f"CGenFF data not found in {cand}: run scripts/fetch_data.sh or set ${TOPPAR_ENV}")
    return cand

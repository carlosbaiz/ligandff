"""LigandFF: CGenFF-compatible atom typing and parameter/charge assignment.

Clean-room reimplementation of the published CGenFF algorithms (docs/SPEC.md). Not the official CGenFF program, and
not affiliated with or endorsed by its developers.
"""
__version__ = "0.1.0.dev0"

TOOL_NAME = "LigandFF"
TARGET_CGENFF_VERSION = "5.0"
DISCLAIMER = "CGenFF-compatible; not the official CGenFF program"


def banner():
    """One-line provenance statement that every output file must carry (SPEC 2.2)."""
    return f"{TOOL_NAME} {__version__} ({DISCLAIMER}); for use with CGenFF version {TARGET_CGENFF_VERSION}"

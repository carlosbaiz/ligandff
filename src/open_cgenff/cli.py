"""Command line: CGenFF-compatible parameters for one molecule (not the official CGenFF program).

usage: open-cgenff INPUT [-o OUT.str] [--json OUT.json] [--name RESI] [--psf OUT.psf] [--gromacs PREFIX]
                         [--toppar DIR]
INPUT: a .mol2/.sdf/.mol file with all hydrogens and bond orders, or a SMILES string (hydrogens added).
Without -o the stream file is printed to standard output. Exit code 1 if typing reported errors.
"""
import argparse
import sys

from . import banner
from .io.stream import to_str, write_report, write_str
from .pipeline import parameterize, read_molecule


def main(argv=None):
    ap = argparse.ArgumentParser(prog="open-cgenff", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="mol2/sdf/mol file or SMILES string")
    ap.add_argument("-o", "--output", help="CHARMM stream file to write (default: standard output)")
    ap.add_argument("--json", help="also write a JSON report (types, charges, penalties, parameter provenance)")
    ap.add_argument("--name", help="residue name (default: molecule name, else LIG)")
    ap.add_argument("--psf", help="also write a CHARMM PSF (for OpenMM's CharmmPsfFile)")
    ap.add_argument("--gromacs", metavar="PREFIX", help="also write PREFIX.top/.gro (needs OpenMM and ParmEd)")
    ap.add_argument("--toppar", help="CHARMM toppar directory (default: data/ from scripts/fetch_data.sh)")
    ap.add_argument("--version", action="version", version=banner())
    args = ap.parse_args(argv)

    result = parameterize(read_molecule(args.input, args.name), args.name, args.toppar)
    if args.output:
        write_str(result, args.output)
    else:
        sys.stdout.write(to_str(result))
    if args.json:
        write_report(result, args.json)
    if args.psf:
        from .io.psf import write_psf
        from .params import default_engine
        write_psf(result, default_engine(args.toppar).db.file.masses, args.psf)
    if args.gromacs:
        from .export import write_gromacs
        write_gromacs(result, args.gromacs, args.toppar)
    for w in result.warnings:
        print(f"warning: {w}", file=sys.stderr)
    for e in result.errors:
        print(f"error: {e}", file=sys.stderr)
    print(f"{result.name}: param penalty {result.param_penalty:.1f}, charge penalty {result.charge_penalty:.1f}, "
          f"{len(result.analogs)} parameters by analogy", file=sys.stderr)
    return 1 if result.errors else 0


if __name__ == "__main__":
    sys.exit(main())

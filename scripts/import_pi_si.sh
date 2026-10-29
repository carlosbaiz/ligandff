#!/usr/bin/env bash
# Unpack the P-I Supporting Information archive (126 "special case" mol2 files; Vanommeslaeghe & MacKerell,
# JCIM 2012, ci300363c_si_001.zip) into data/pi_si/ (git-ignored; copyrighted SI is never committed, SPEC 2.3).
# The archive is not fetched automatically: download it from the journal and pass its path.
# usage: bash scripts/import_pi_si.sh /path/to/ci300363c_si_001.zip
set -euo pipefail
cd "$(dirname "$0")/.."
[ $# -eq 1 ] || { echo "usage: $0 ci300363c_si_001.zip" >&2; exit 2; }
mkdir -p data/pi_si
python -I -c "import sys, zipfile; z = zipfile.ZipFile(sys.argv[1]); z.extractall(sys.argv[2]); print(len(z.namelist()), 'files')" "$1" data/pi_si
ls data/pi_si/specialcases/*.mol2 | wc -l

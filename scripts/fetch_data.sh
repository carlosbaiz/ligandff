#!/usr/bin/env bash
# Fetch the public CGenFF data at pinned commits into data/ (not committed; see docs/SPEC.md section 2.3).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data

fetch() {  # url dir commit
  if [ ! -d "data/$2/.git" ]; then git clone --filter=blob:none "$1" "data/$2"; fi
  git -C "data/$2" fetch -q origin "$3" 2>/dev/null || true
  git -C "data/$2" checkout -q "$3"
  echo "data/$2 @ $(git -C "data/$2" rev-parse --short HEAD)"
}

# CHARMM36 + CGenFF 5.0 force-field files (Feb 2026 release): top_all36_cgenff.rtf, par_all36_cgenff.prm,
# par_all36_cgenff_increments.prm live in charmm/toppar/toppar/
fetch https://github.com/Lemkul-Lab/toppar_c36_releases.git toppar_c36_releases ed3bbdaa4eb0a0152e577480bb0a9d22cd9743c3
# CGenFF 5.0 paper data: 1,390 training molecules (mol2, SMILES), optimized parameters, validation data
fetch https://github.com/mackerell-lab/cgenff5.git cgenff5 f3e28e93a9aa12a9c8210e6dfd91808605213f29

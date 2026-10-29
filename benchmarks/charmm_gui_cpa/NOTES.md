# B3: CHARMM-GUI topologies (cryoprotectants and small ions)

Source: Carlos Baiz, 2026-10-29 (`SOURCE_README.txt`). Generated with CHARMM-GUI v1.7 (`psf2itp_mol.py`).
`raw/<A>.itp` = molecule topology (atoms with CGenFF types and charges, bonds; no bond orders);
`raw/cgen_<A>.itp` = parameters (`[ atomtypes ]`, bonded types).

Known issues (do not "fix" the files; handle in code):
- `glycerol`: CHARMM36 carbohydrate types (CC322, HCA2, OC311 …), not CGenFF → excluded from CGenFF scoring.
- `mbn`: `mbn.itp` and `cgen_mbn.itp` contents are swapped. Charges are custom ("charges from Min Cho's paper") and sum
  to 0 although S is a thiolate (SG302) → typing scored with net charge −1; excluded from charge comparison.
- `tip3p`, `tip4pew`, `tip4pice`, `au`: water models and gold, not CGenFF ligands → excluded.
- Version: parameters match CGenFF 4.6 wherever 4.6 has the term ⇒ charges come from a 4.x-era program whose increments
  are not public. Use charges only as a version-drift reference, not as an exactness test.
- Duplicates: `prgl`, `proh`, `otpr` (propylene glycol) and `buta`, `ttbd` (2,3-butanediol) are the same molecules.

Scored species: 34; current pilot result: 34/34 molecules, 545/545 atoms (FINDINGS §3.4).

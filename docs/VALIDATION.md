# Validation report (M7)

LigandFF 0.1.0, a CGenFF-compatible reimplementation (not the official CGenFF program), targeting CGenFF 5.0.
Date: 2026-10-29. Reproduce with `python -m ligandff.benchmarks.validation` (writes `build/report/validation.md`,
`report.md`, `report.json`; CI runs it on every push). Benchmarks and splits: SPEC section 7. Test splits are
reported as aggregates only (SPEC 7.2).

## 1. Definition of done (SPEC section 1)

| target | measured | met |
|---|---|---|
| types B1-test: molecules >= 96 %, atoms >= 99.5 % | 97.35 % (220/226), 99.69 % | yes |
| types B2-test: molecules >= 95 %, atoms >= 99.5 % | 97.80 % (178/182), 99.62 % | yes |
| types B3: 100 % | 100 % (34/34), 545/545 atoms | yes |
| charges B0: exact (charges and penalties to 0.001 e) | 3/3 | yes |
| charges QT unchanged (2,255 residues, median RMSE 0.0098 e) | 2,255 residues, 0.0098 e | yes |
| formal charges QF >= 887/218/669/166 | 887/218/669/166 | yes |
| impropers QI >= 884/219/705/175 | 884/219/705/175 | yes |
| analogy B4-test: same source >= 85 %, TPS within ±0.5 >= 80 % | 92.1 %, 92.5 % (1,447 cases) | yes |
| output: `.str` + PSF build an OpenMM system for every B1/B2/B3 molecule | 2,074/2,074 | yes |

All regression floors of SPEC 7.4 are met.

## 2. Train and test side by side

| benchmark | train / dev | test |
|---|---|---|
| B1 types, molecules fully correct | 99.00 % (895/904) | 97.35 % (220/226) |
| B1 types, atoms | 99.84 % | 99.69 % |
| B2 types, molecules fully correct | 98.90 % (720/728) | 97.80 % (178/182) |
| B2 types, atoms | 99.74 % | 99.62 % |
| QF formal charges (molecules with correct types and complete increments) | B1 887/888, B2 669/670 | B1 218/218, B2 166/166 |
| QI improper centres (molecules with correct types) | B1 884/895, B2 705/720 | B1 219/220, B2 175/178 |
| B4 same source / TPS within ±0.5 | 93.0 % / 92.6 % (5,787) | 92.1 % / 92.5 % (1,447) |

The train-test gap is 1–1.7 points for typing and under 1 point for the analogy: no sign of strong overfitting, but B1
and B2 train and test come from the same population (CGenFF training residues).

## 3. End-to-end output (set OUT)

Every B1/B2/B3 molecule through typing, bonded parameters, charges and lone pairs, written as `.str` + PSF and read
by OpenMM with the CGenFF files: one particle per atom and lone pair, finite energy at an RDKit starting geometry.

| set/split | molecules | build OK | lone pairs | need analogs | param penalty < 10 | 10–50 | > 50 |
|---|---|---|---|---|---|---|---|
| B1/train | 904 | 904 | 0 | 15 | 894 | 0 | 10 |
| B1/test | 226 | 226 | 0 | 5 | 222 | 0 | 4 |
| B2/train | 728 | 728 | 24 | 10 | 720 | 0 | 8 |
| B2/test | 182 | 182 | 7 | 3 | 179 | 0 | 3 |
| B3 | 34 | 34 | 0 | 14 | 33 | 0 | 1 |

Throughput on one core (paper/collect.py, force-field files loaded once): median 2 ms, mean 6.4 ms per molecule for the
full pipeline (~160 molecules/s); molecules that need analog searches take 10–900 ms. Typing alone: 200–370
molecules/s. The OpenMM build adds ~0.5 s per molecule (reading the CGenFF files each time). The OUT table's rate is
measured inside the parallel OpenMM workers and is lower because of contention. The P-II figure for
the official program is ~30 molecules/s on one 2012 core.

Found and fixed during M7 (FINDINGS §11):
- 265 molecules failed to build because starting coordinates were embedded from connectivity alone (an N+ or a
  phosphate P with only single bonds fails RDKit's valence check). Coordinates now come from the input molecule with
  its bond orders and formal charges.
- Dihedrals through a linear centre (C≡N, C≡C) received analogs with penalties up to ~3,000 in correctly typed
  training molecules (acetonitrile, the alkynes). CGenFF defines only two such dihedrals and the program recorded no
  analog for any, so they are now omitted when no exact parameter exists. Molecules needing analogs: 68 → 47; with a
  parameter penalty above 50: 47 → 26 (on the development splits all but two of these, PHOS and 2300YC, are molecules
  whose types we get wrong).

## 4. Limits

1. **No measurement on new molecules (D4).** Every reference is a CGenFF training molecule (B1, B2), a CHARMM-GUI
   topology of the 4.x era (B3), the program's own analogies stored in the parameter file (B4), or the three molecules
   printed in the P-v5 SI (B0). How closely the tool matches the official program on new drug-like molecules is not
   measured. B1/B2 bound typing accuracy from above: the rules were written on these molecules (train splits).
2. **Charges of molecules that need analogs** are validated only on B0 (M064, one analog). QT checks the increment
   model on complete training residues; the charge-analogy path (nonbonded matrix, dihedral cap at 50) is not
   validated beyond B0.
3. **Penalty values.** The fitted matrix (M5) ranks candidates well (B4 same source 92 %), but substitutions never
   observed in the parameter file cost 20 × the relation median. Penalties built from them (typically > 100, up to
   thousands) are a "do not trust" signal, not values the official program would print. Penalties below ~10 come
   from observed elements and agree with the program on B4 (TPS within ±0.5 for 92 %).
4. **Linear-centre dihedrals** are left out of our PSF and OpenMM system (section 3). CHARMM itself may still generate
   them from the residue (`AUTOGENERATE ... DIHE`), exactly as for official training residues such as ACN; how CHARMM
   treats the missing parameters there is not checked. The charge side is unchanged: dihedral increments through a
   linear centre come from analogs above the cap of 50, so they add no charge but a small charge penalty
   (acetonitrile: 4.3).
5. **Typing residuals** (dev splits): mostly zwitterionic/mesoionic resonance structures (e.g. 2300VB, PNTM),
   exotic heterocycles and alternation `altnum` choices; listed in the benchmark report's confusion tables.
6. **Order dependence:** resonance ties and `altnum` can depend on input atom order (SPEC 5.3, 5.5).
7. **Formats:** GROMACS export refuses halogen lone pairs (ParmEd does not convert them; use `cgenff_charmm2gmx`).
   mol2 input carries no formal charges; they are derived from valences.
8. **Version drift:** B3 charges come from a 4.x-era program and are not used as an exactness test; B4 annotations
   mix program versions, part of the remaining 8 % misses.

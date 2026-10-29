# Findings so far (feasibility study and pilot, 2026-10-28/29)

Everything here was measured with the code in `prototype/` on the data pinned in `scripts/fetch_data.sh`.
Reproduce with the commands in `prototype/README.md`.

## 1. Feasibility in one paragraph
The CGenFF program's algorithms are fully published (P-I, P-II, P-X), and the data it trains on (topology, parameters,
and, since CGenFF 5.0, the charge increments) are public. Only two hand-built data sets are unpublished: the atom-typing
rule file and the substitution-penalty rule file. The charge half of the program has been reproduced exactly; a
pilot typing rule set reached 93.8 % of held-out molecules fully correct after three refinement rounds. The remaining
effort is concentrated in typing rules (long tail) and the penalty matrix.

## 2. Charges

### 2.1 Exact reproduction of program 5.0 output
`prototype/charges.py fixtures` on the three molecules printed with program-5.0 output in the P-v5 SI (Table S10):

| molecule | needs analogy | charges | per-atom charge penalties |
|---|---|---|---|
| M064 (diisopropyl ketone) | 1 dihedral, CG311-CG2O5-CG311-CG331 from CG321-CG2O5-CG311-CG331, penalty 0.6 | exact (0.001 e) | exact (0.334 / 0.657 / 0.669 …) |
| M166 (cyclopentanone) | none | exact | exact (all 0) |
| M223 (indane) | none | exact | exact (all 0) |

Two implementation details were settled by this test: (a) for the missing dihedral the program took the charge increments
from the same analog parameter as the bonded terms (here the substitution is at an outer position, so the bonded and
nonbonded matrices agree; the general rule may differ); (b) in P-II eq 4, β_n is the net charge one parameter instance
moves onto the atom.

### 2.2 Formula check on all training residues
For 2,255 RTF residues (neutral and charged, no lone pairs or inter-residue bonds) whose increments are all present,
increment charges + inferred integer formal charges reproduce the hand-optimized topology charges with median
per-molecule RMSE 0.0098 e; 76.0 % of atoms within 0.01 e, 94.0 % within 0.05 e (fit residual of the increment model;
P-II reported 0.0082 e RMSD for the 2012 fit). Inferred formal charges match P-I's description: −1 on the carboxylate
carbon (CG2O3, 133/135), +1 on amidinium/guanidinium carbons (CG2N2, CG2N1), +1 on ammonium N (NG3P0-3), +1 on some
NG2R61 (pyridinium-type) and CG2R53 (imidazolium C2), −1 PG1, −2 PG2, −1 OG312/SG302/SG3O1.
The full list of inferred (type, q0) pairs is in the benchmark report (set QT, M1). Other frequent pairs not covered
above: CG2R63 +1 (36 atoms), OG2D1 −1 (19), OG2D2 −1 (12), SG2D1 −1 (11), CG2R64 +1 (10); plus a tail of 1–7-atom
pairs (e.g. CG2R61 ±1) that may be fit residuals > 0.5 e rather than rule charges. To be resolved in M3 (SPEC 5.7).

### 2.3 Halogen lone pairs
P-X: LP charge +0.05 e fixed, no LJ, colinear, 1.64/1.89/2.20 Å for Cl/Br/I (matches the RTF). Halogen charge offset is
not published; on 38 halogen-LP pairs, "halogen = increment charge − 0.05" gives max error 0.04 e (mean 0.017 e);
"halogen = increment charge" gives max 0.06 e (mean 0.035 e).

## 3. Atom typing pilot

### 3.1 Benchmark B1
RTF residues with SMILES comments, mapped atom-by-atom onto RDKit molecules by graph isomorphism (only exact bond-graph
matches kept): 1,130 molecules, 27,512 atoms, 126 types; all are v5.0 ZINC-derived training compounds (the 937 older model
compounds have no SMILES; see B2 in SPEC §7). 94 SMILES comments had lost bond orders (radical atoms) and were repaired
from connectivity + net charge (1 unrepairable); 264 residues were skipped (atom-count or bond-graph mismatch).
Split 80/20 with seed 42. A naive lookup on each atom's immediate environment types ~89 % of atoms but only ~18 % of
molecules fully correctly: whole-molecule fidelity needs explicit rules.

### 3.2 Learning curve (only training-split failures were inspected while writing v1–v4)

| rule set | train molecules fully correct | held-out molecules | held-out atoms |
|---|---|---|---|
| v1: written from MASS-line descriptions + P-I only | 50.1 % | 50.4 % | 94.9 % |
| v2: + 1 round of failure analysis | 85.3 % | 88.5 % | 99.0 % |
| v3: + 2 rounds | 92.9 % | 92.9 % | 99.4 % |
| v4: + 3 rounds (= `prototype/typer_pilot.py`, ~110 rules) | 96.7 % | 93.8 % | 99.2 % |

Wall time for building the benchmark and v1–v4: ~13 minutes of agent time. After v4, held-out failures were inspected
once to characterize them (so B1-test is now semi-sealed): ring-embedded amidinium/guanidinium cations, an amine N-oxide
(neighbouring CH2/CH typed as "adjacent to positive N"), a sulfonic acid (SG3O2, not SG3O1), a cyclic sulfinamide, and
5/7-fused tropolone-like systems (CG2RC7).

### 3.3 CGenFF conventions learned from data (not obvious from type descriptions)
- CG2RC0: carbon of an aromatic 6-ring that is also in any 5-ring (incl. non-aromatic, e.g. indane) or 4-ring
  (P-v5 rule); NOT for 5/5-fused systems (those keep CG2R5x types). NG2RC0 needs the 5-ring to be planar.
- CG2R67 (biphenyl) also when the aryl–aryl bond lies in a non-aromatic ring (fluorene-, carbazole-, phenanthridine-like).
- CG2R57: 5-ring sp2 carbon single-bonded to an atom whose smallest ring is 5-membered (aromatic or not), not only bipyrroles.
- CG2R51/52/53 for C=C/C=N inside any 5-ring, incl. non-aromatic (3-pyrroline, maleimide); P-I's `sub NG2R5: or (arom 5) (ring2 5) (ring23 5)` agrees.
- CG2R53 also for carbonyl carbons in 5-ring lactones/lactams and 4-ring lactams; cyclobutanone carbonyl is CG2O4 (with OG2D3).
- CG2R64: aromatic 6-ring carbon with ≥ 2 N neighbours (ring or exocyclic).
- CG322/CG312/CG302 by number of F (1/2/3) regardless of H count; 5-ring CH with F stays CG3C51.
- CG3RC1 only for true fusion/bridgehead atoms of a ring system containing a ring ≤ 5; spiro atoms are not CG3RC1.
- NG3C51: any non-amide sp3-type N in a 5-ring (incl. tertiary, N-aryl, isoxazolidine N), checked before hydrazine NG3N1.
- NG2S3 only for NH2 on an aromatic 6-ring atom; aryl-NH-R / aryl-NR2 are NG311/NG301; NH2 on aromatic 5-ring C → NG321;
  NH2 on cationic aromatic N → NG3N1.
- NG2R53 only when N and the carbonyl C share a 5-ring; NG2R43 for 4-ring lactams.
- OG3C51 also for oxetane O.
- H: HGR62 on 6-ring aromatic C adjacent to N or to the C–F carbon (CG2R66), and on CG2R62/CG2R64; HGR63 when the ring is
  cationic; HGP5 on C attached to NG3P0; HGPAM1/2 on NG311/NG321 unless N is bonded to S/O/N (then HGP1); HGP1 on NG3C51;
  HGP4 on NG2S3; HGP2 on N of cationic systems (ammonium, amidinium, imidazolium).
- "Positive N" for CG314/CG324/CG334/CG3C53/CG3C54 includes aromatic N of cationic rings and amidinium N.
- Resonance tie-breaking matters (quinolinium: charge on ring N vs exocyclic N have equal P-I scores).

### 3.4 Benchmark B3 (CHARMM-GUI, supplied by Carlos)
39 files → 34 CGenFF-typed species scored (excluded: 3 water models, gold, glycerol with CHARMM36 carbohydrate types).
26 of the 34 are not CGenFF reference residues. Pilot v4: 545/545 atoms, 34/34 molecules; v1 (documentation only):
33/34 (guanidinium wrong). The set is chemically simple (polyols, glycol ethers, amides, small ions). Parameters in these
files match CGenFF 4.6 wherever 4.6 has the term (1,403/1,403 dihedrals, 276/276 bonds) ⇒ generated with a 4.x-era
program; its increments are not public. For the 20 molecules fully covered by exact 5.0 increments, 5.0 charges equal
the CHARMM-GUI charges for 2 (NMA, pinacol), are within 0.01 e for 11, max 0.069 e: version drift, not a formula error.
MBN: `mbn.itp`/`cgen_mbn.itp` contents swapped; charges custom (net 0 despite a thiolate S): typing only.

## 4. Penalty matrix evidence
`prototype/penalty_annotations.py` on `par_all36_cgenff_increments.prm`: 7,234 distinct annotated analog parameters;
5,089 differ from their source at one position (156 bonds, 1,518 angles, 3,302 dihedrals, 113 impropers) → 107
bonded-matrix and 280 nonbonded-matrix elements over 108 types (≈ 1 % of each ~27,000-element matrix); 106/107 and
279/280 are internally consistent. P-II Table 1 matches the bonded matrix exactly (NG311→NG301 0.5, NG311→NG321 1.0,
NG301→NG311 0.5, NG3C51→NG311 2.0); the nonbonded matrix differs (e.g. NG3C51→NG311 3.0, NG321→NG311 1.2).

## 5. Compute
All benchmarks run in seconds on one CPU core (RDKit). Projected cost of completing typing + penalty matrix with an
agent loop: ~150–400 refinement rounds, ~30–160 M input tokens (mostly cache reads) + 1–4 M output tokens.

## 6. Perception (M2)
Production perception (`src/open_cgenff/perception/`) on the P-I material and the pilot typer:
- P-I worked examples reproduced exactly: Fig. 4 yields structures H, I, J with scores 12, 11, 11 in that order and
  keeps I; pyridinium yields two structures scoring 6 with the pre-0.9.7 coefficients (as printed in case study 1)
  and 11 with the current ones. The choice order (lowest-index neighbour, double before single, neutral before
  charged) is read from the Fig. 4 trails.
- All 126 P-I SI special cases run (1.9 s total). Two warn: `fuller_all_ar` (search budget reached; P-I notes
  fullerenes are expensive) and `isoindenate` (no valid structure without a carbanion).
- Score scope: counting the whole molecule's charge in `Z_mol` made B2-train worse (621 vs 625 molecules, −93 atoms):
  a fixed phosphate dianion makes an N+ that breaks purine aromaticity cheaper (AMP, GTP analogues). Counting only
  the charges chosen on the subgraph gives 627/728 and equal B2-test aggregates. Adopted as the interpretation.
- Pilot typer on production vs its own perception: B1-train 873 vs 874, B1-test 213 vs 212, B2-train 627 vs 625,
  B2-test 158 vs 158, B3 34/34; atoms correct B1-train 21,999 vs 21,995. B1-train changes: 2300RV and 230051 fixed;
  2300J3 and 2300PE lost because their radical-repaired SMILES put the enolate charge on the other oxygen and P-I
  takes non-aromatic input bond orders literally (the pilot re-chose among RDKit resonance forms); 230003 lost because
  a pilot rule relied on RDKit moving the amidinium charge through non-aromatic conjugation. Remaining B2-train
  differences come from B2 inputs whose RDKit bond orders are unusual (PNTM gets an `[O+]=` oxonium, SM033 an azide
  form): see `benchmarks/b2_rtf_legacy.NOTES.md`.

## 7. Typing rules (M3)
Production typer = P-I rule language (`src/open_cgenff/typing/`) + rule file v1 (`typing/rules/cgenff_v5.rules`).
Only B1-train/B2-train failures were inspected; test splits are aggregates.
- Port: the pilot rewritten in the rule language reproduces the pilot on production perception atom for atom on
  B1/B2-train (873/627), after one interpretation was fixed: nested `ne` must not revisit the atom it came from.
  P-I's printed `NG2R62` rule (`ne ( inring arom 6 or ( el N ) ( ne ( inring el N arom 6 ) ) )`) only works that way,
  so predicates written "from the other side" (amidinium C seen from its N; nitro O seeing its N) were rewritten.
- P-I's `NG2R62` (another ring N ortho/meta in the same aromatic 6-ring) beats the pilot's broader rule (any aromatic
  heteroatom within two bonds): +2 B1-train, +3 B1-test, +1 B2-train, nothing lost.
- Conventions learned on B2-train (counts = train molecules/atoms with the pattern):
  H on an aromatic C next to a C–halogen is HGR62 (50/50 atoms, 24 halobenzenes); sulfate/sulfonate S is SG3O1 only
  with an O− (neutral sulfonic acids and esters SG3O2); thio/dithio-phosphate S is SG2P1/SG2P2 and the terminal O on
  such a P is OG2S1; N–P and N–SO3− are NG2S3 with HGP1; methylamine/dimethylamine/trimethylamine methyls are
  CG3AM2/1/0 with HGAAM*; CH3 on any N+ incl. N-oxides and nitro is CG334; N-oxide O− is OG312 (nitro O stays OG2N1);
  boronate O− is OG2D2; imidazolium C2-H is HGR53 (also in imidazolinium); H on C=N+ and on formate is HGR52;
  N-aryl on an aromatic 6-ring N is CG2R67 like biaryls; azide N1/N2 are NG1T1; isocyanate C is CG2O7; carbonate C
  CG2O6; beta-lactone carbonyl CG2O4 (4-ring CG2R53 needs N); conjugation for CG2DC* counts C=X but not C≡X; an
  amidinium inside an aromatic 6-ring (protonated cytosines: ring N next to a C with an exocyclic N–H, system
  cationic) gives NG2P1 ring/exocyclic N, CG2R64 carbons and a ring that is *not* cationic for CG2R62/HGR63;
  pyrimidinium C2 (cationic ring, between two N) is CG2R62, not CG2R64; imidates `N=C–O−` are typed as amides
  (CG2O1/OG2D1); C=Se is treated like C=O (CG2R63).
- B1-train: 7-ring sp2 carbons are CG2R71 (incl. C=O) and 7/5 fusion carbons CG2RC7 only when the 7-ring is aromatic
  and fused to a planar 5-ring (azulene-like; benzodiazepine-like 5/7 fusions keep CG2R5x); a 5-ring N bonded to an
  aromatic 5-ring is NG2R57 whatever its own ring class (lactam N, indoline N); 2-fluoropyrimidine C is CG2R64 (N wins
  over F); exocyclic C=N on a 5-ring C is CG25C1 (the double-bond partner does not count as a heteroatom); isourea and
  isothiourea C are CG2N2.
- Remaining B1-train failures (9): 2300J3, 2300PE, 23008L, 2300VB, 230073 have input structures from the radical
  repair or zwitterionic SMILES that P-I takes literally (FINDINGS §6); 23001K, 23002G, 2300RO are exocyclic-double-
  bond conventions with one example each and conflicting neighbours; 2300CY (7-ring fused to pyrazole typed
  azulene-like although not aromatic by P-I counting). B2-train (8): PNTM, B1MA, PEGM, SM212 input artifacts; small
  S–O/S–S molecules SM006, SM053, SM055, SM152.
- Formal charges (SPEC 5.7, set QF): with the rule charges, 887/888 B1-train, 218/218 B1-test, 669/670 B2-train and
  166/166 B2-test molecules (complete increments, no lone pairs, types correct) have exactly round(q_ref − q_incr).
  Exceptions: 23012O (C=O carbon −1, O +1: a dipole the increments do not capture) and C34H.
- B2 input fix: DetermineBondOrders produced charge-separated structures for 14 residues (thioesters as `[S+]=C–[O−]`,
  2-alkylthio-pyrimidinones, a thiosulfate, a thiophosphate). `b2.py` now flips alternating paths between a +1 atom
  and a −1 atom. Pilot baseline on the rebuilt B2: 633/728 train, 161/182 test (was 625/158).

## 8. Parameters by analogy (M4)
`src/open_cgenff/params/` on `par_all36_cgenff.prm` (v5.0): 1,305 bonds, 6,776 angles (721 with Urey–Bradley),
19,629 dihedrals (29,044 lines; up to 6 terms), 489 impropers; 7,234 `from …, penalty=` annotations.
- Ordering: the header's rules (element priority C < N < O < P < S < halogens < misc < H, alphabetical within an
  element; section sort keys bonds 1,2 / angles 2,1,3 / dihedrals 2,3,1,4 / impropers 1,4,2,3) reproduce the file:
  every bond, angle and dihedral is in canonical orientation, and the line order is the collated order except for two
  misplaced lines (dihedral NG2D1 CG25C1 CG25C2 OG2D3, improper NG2S3 CG2R64 HGP4 HGP4).
- M064 (SI Table S10): with the provisional matrix the search returns the program's analog (CG321-CG2O5-CG311-CG331,
  TPS 0.6: one outer substitution CG311→CG321 = nonbonded element 0.6, observed in the annotations); with the
  nonbonded-everywhere analog for increments, charges and per-atom charge penalties match the official output
  (charge penalty 0.669). M166/M223: no analogs, both penalties 0.
- Provisional matrix (SPEC 5.10 status): 106 bonded / 260 nonbonded observed elements after excluding annotations with
  a bond-group difference. Fallback medians (bonded / nonbonded): same family with 2/1/0 further common characters
  0.6/1.0/1.5 and 2.0/3.0/4.0; same element, other family 9.0 / 68.0 (H types: 0.2 bonded, from HGP5→HGA1/HGA2);
  different elements 9.0 (unobserved: maximum observed element) / 46.5. The nonbonded matrix makes element changes
  cheaper than CG2↔CG3 changes (CG2N2→CG321 69, NG2R52→CG2R51 46.5): data, not a bug.
- Leave-one-out baseline (every 10th annotated parameter, 725; the matrix elements come from the same annotations, so
  this is optimistic): same source 36.4 % (bonds 82 %, angles 33 %, dihedrals 37 %, impropers 33 %), |ΔTPS| ≤ 0.5
  46.3 %. Most misses involve fallback elements or multi-substitution analogs: the M5 problem.
- Coverage: with reference types, every bond and angle of the 904 B1-train residues exists in the file; all dihedrals
  except 22 through a linear centre (X–CG1T1≡CG1T2–Y, X–CG1N1≡NG1T1 neighbourhoods, 7 molecules) that CGenFF does
  not define. Of 632 improper instances from the rules, 619 match a parameter exactly (outer atoms in any order).
- Speed: 0.02–0.2 s per drug-like molecule for all bonded analogs plus charge analogs (search indexed by inner types).
- Impropers (SPEC 5.13): the rule file v1 carried no `impr`. Ground truth = IMPR lines of the B1/B2 residues (first atom
  always central, 617/617 on B1-train). Centres follow mostly from chemistry, not type: every exocyclic C=O/C=S carbon
  (CG2O*, CG2R63, CG2R53 with C=O: 100 %), CG2N1/CG2N2, sp2 C with an exocyclic C=N (CG2D1, CG2DC*, CG25C1,
  CG2R61, CG2R64), CG2D1O with an N neighbour (but ring C=C–O carbons have none, against P-I's N/O/S condition),
  ring C next to a ring N carrying an exocyclic N (CG2R64 90/100, CG2R53 4/4), aryl NH2 N (NG2S3, 96/98).
  Remaining B1-train mismatches: CG2R64 with an exocyclic N but no IMPR in the residue (8; e.g. 23001G).

## 9. Output and export (M6)
- Halogen lone pairs: the v5.0 RTF has 66 LONEPAIR COLINEAR lines, all at 1.64/1.89/2.20 Å (Cl/Br/I; no DIST 0.0),
  two of them on aliphatic halogens (SM200 BRGA1, SM211 CLGA1; not generated: P-X puts lone pairs on aromatic
  halogens only). On the B2-train residues with aromatic-halogen lone pairs, halogen charge − increment charge is −0.05 e
  for most monohalides (CHLB/PCLT/PMHA −0.06, BROB −0.049, IODB/SM196–SM205 −0.05) and between −0.04 and +0.03 e in
  polyhalobenzenes, whose other atoms are also hand-adjusted (up to 0.1 e): "halogen = increment − 0.05" stays.
  Hygiene note: the calibration script also printed the per-residue rows of 7 B2-test residues (dihalo-/trihalobenzenes
  and SM207); they were not used for any decision, and the rule was not changed.
- OpenMM: our PSF (XPLOR types, NUMLP section) + `CharmmParameterSet(rtf, prm, .str)` builds systems for all 34 B3
  molecules and all 31 B1/B2 molecules with lone pairs; the virtual site lands 1.64 Å beyond Cl on the C–Cl axis.
  The CGenFF parameter file defines LPH with zero LJ. Pipeline on all 2,040 B1/B2 molecules: no crash, no typing
  error; 54 molecules need at least one analog parameter (training residues whose types we get wrong, or the
  linear-centre dihedrals of FINDINGS §8). ~0.5 s per molecule including the OpenMM build.

# Rules changelog

One entry per typing- or penalty-rule change: date, rule(s), motivating molecules (dev splits only), metric deltas
(B1-train/test, B2-train/test, B3; molecules fully correct).

## Pilot (prototype/typer_pilot.py; Python rules, not the production rule language)
- v1 (documentation only): B1 train 50.1 %, test 50.4 %.
- v2: CG2RC0 for 6/5 fusions incl. non-aromatic 5-rings; spiro vs fused bridgeheads; resonance-aware cationic N;
  5-ring sp2 types in non-aromatic rings; NG2S3 restricted to Ar-NH2; H types by parent. → 85.3 % / 88.5 %.
- v3: CG2R57 generalized; CG2R67 for in-ring biaryl links; NH2 on 5-ring atoms; cationic definitions from best
  resonance structure; cyclobutanone; oxetane O; vinyl F; conjugation detected from either end of C=C. → 92.9 % / 92.9 %.
- v4: biaryl test restricted to shared aromatic 6-rings; NG2R62 via aromatic heteroatoms; NG2RC0 needs planar 5-ring;
  F counts after 5-ring check; resonance tie-break prefers ring-centred cation; bridgehead definition. → 96.7 % / 93.8 %.
Details of each convention: `docs/FINDINGS.md` §3.3.

## Perception (M2, no rule changes)
- Pilot typer on production perception (`benchmarks/plugins/pilot_on_production_perception.py`): B1 train 873/904
  (−1), test 213/226 (+1); B2 train 627/728 (+2), test 158/182 (=); B3 34/34. Explanations: FINDINGS §6.
- B2 baseline (pilot typer, own perception): train 625/728 (85.85 %), test 158/182 (86.81 %); now a CI floor.

## M3: production rule language and rule file v1 (2026-10-29)
Counts are molecules fully correct, B1-train/B1-test/B2-train/B2-test/B3; test splits are aggregates only.
Before: pilot on production perception 873/213/627/158/34 (M2 B2 build).

1. Port of the pilot to `cgenff_v5.rules` (language: SPEC 5.5). Fixes during the port: amidinium predicate written as
   seen from its N (nested `ne` does not revisit the parent; 23001V, 2300DS, BAMI, 2IMP); `fusion` keyword with the
   pilot's definition (CG3RC1; 23005K); NH4+ is NG3P3; CG2N2 for an amidine conjugated with an amidinium N (2300CU).
   → 873/213/627/158/34 (identical to the pilot on production perception).
2. HGR62 next to an aromatic C–halogen (CHLB, 12DCB, … 24 molecules); SG3O1 needs an O− (MMST, MSMT, PSMT, SM149);
   SG2P1/SG2P2/OG2S1 thiophosphates (GTPG, GTSS, GTGN); NG2S3 + HGP1 on N–P (PHA, TPHC, THNP, SM150, SM156);
   CG3AM*/HGAAM* (MAM1, DMAM); CG334 on N-oxide/nitro N+ (TMAO, NIME); OG312 for N-oxide O− (TMAO); OG2D2 for
   boronate O− (BONN, BONE, …); HGR53 imidazolium C2-H (2IMP, B7MG, …); HGR52 on C=N+ and formate (SCH2, SCH3, FORA);
   CG2R67 for N-aryl (1PH2PO, 3PHURA, 43PPY); NG1T1 azide N1 (SM033, SM217); FGP1/ALG1 (ALF4); CG323 for CH3–S−
   (MES1); SG301 for S–O (SM006, SM052, SM055). Nitro O kept OG2N1 (NITB). → 879/215/691/173/34.
3. B2 builder: charge-separation fix (`b2.py:_neutralize`, 14 residues; not a rule change). → B2 691→698 train,
   173→176 test with unchanged rules; pilot baseline on the new build 633/161.
4. Protonated cytosine-type ring amidinium: NG2P1 ring and exocyclic N, ring not cationic (B3MC, NCYP, PNCP, BK2C,
   BR2C); pyrimidinium C2 CG2R62 (PYRH, PYRF); imidate anion as amide (ABMB, ABEB, ABNB); isothiourea CG2N2 (DH3T,
   23002W); C=Se like C=O (BSEU); carbonate CG2O6 + OG2D2 (CO3); isocyanate CG2O7 (MICY); beta-lactone CG2O4/OG2D1
   (2OXT); no conjugation through C≡C (BEYN); azide N2; sulfamate NG2S3 (NMSM, NESM, 2300NR); N–S H is HGP1.
   → 881/217/718/178/34.
5. Azulene-type 7/5 systems: CG2RC7 fusion C and CG2R71 (incl. C=O, OG2D3) when the 7-ring is aromatic and fused to a
   planar 5-ring (23001Q, 2300DN, AZUL); a first version without the aromatic condition broke 230028, 23009I, 23009L,
   2300FW, 2300NR and was discarded. NG2P1 for N+ with an exocyclic C=N+ on a 5-ring, CG334 next to it (2300YF).
   → 884/217/719/178/34.
6. NG2R57 for a 5-ring N on an aromatic 5-ring, any ring class (2300J6, 23012A); CG2R64 for 2-fluoropyrimidine
   (2300L4); CG2R57 link may close a 6-ring (2300P5); exocyclic C=N partner is not a heteroatom for CG25?O (2300YW).
   → 889/217/719/178/34.
7. NG2R62 replaced by P-I's published rule (ring N ortho/meta in the same aromatic 6-ring) instead of the pilot's
   "aromatic heteroatom within two bonds" (230071, 2300YC, PIPA). → 891/220/720/178/34.
8. CG2N2 for isourea (2300CK) and for an amidine bonded to an amidinium C (23003U); exocyclic C=C on a P-I-aromatic
   5-ring is still CG25?O (23006I); conjugation through an amidinium N for CG2DC* (230003). → 895/220/720/178/34.
9. Formal charges (SPEC 5.7; no type changes): P-I's `sub CG2D : ne ( bo 2 el N nb 4 ) charge 1` and the same for
   CG2R52/CG2R53; CG2R64 +1 with an exocyclic C=N+; PG1 −1, PG2 −2, SG3O1 −1, SG302 −1, OG312 −1, boronate OG2D2 −1,
   carbonate CG2O6 −2, ALG1 −1, cationic NG2R67/NG2RC0 +1. QF: 887/218/669/166 molecules with correct formal charges.

After M3: 895/904 (99.0 %), 220/226 (97.3 %), 720/728 (98.9 %), 178/182 (97.8 %), 34/34; these are the new floors.

## M4: impropers (2026-10-29)
10. `impr` (SPEC 5.13; no type changes): `sub CGO … impr` (exocyclic C=O/C=S); pass-through `sub CG2 : … impr` routes
    before the plain `sub CG2`: P-I's exocyclic C=N and C=N + exocyclic N conditions, P-I's third condition narrowed to
    an N neighbour (plus O/S only for acyclic atoms: PY01, 23007W, SM172 have ring C=C–O without IMPR), aromatic C with
    a ring N and an exocyclic N (2-aminopyrimidine type; 2300B7, CSGU); `typ NG2S3 … ne ( arom 6 ) impr` (230001).
    QI (molecules with exactly the reference improper centres, B1-train/test, B2-train/test): 467/108/407/103 →
    884/219/705/175. Typing unchanged: 895/220/720/178/34; QF unchanged.

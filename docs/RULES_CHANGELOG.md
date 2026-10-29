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

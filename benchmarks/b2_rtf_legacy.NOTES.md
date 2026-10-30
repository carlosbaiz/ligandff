# B2: legacy RTF model compounds (no SMILES in the RTF)

Built by `src/ligandff/benchmarks/b2.py` from `top_all36_cgenff.rtf`; only residue names are stored here
(`b2_rtf_legacy.split.json`, 728 train / 182 test, seed 42). B2-test is sealed: aggregate metrics only.

Input quality (bond orders come from RDKit DetermineBondOrders on connectivity + net charge, not from the RTF);
issues seen on B2-train while checking perception in M2 (handle in the builder, never by editing reference types):
- PNTM (pentamidine-like bis-amidinium): one ring gets a quinoid structure with an `[O+]=` oxonium ether.
- SM033 (ethyl azide): azide drawn as `N=[N+]=[N-]` (typed correctly since M3).
- Since M3 the builder removes charge separation along alternating paths (`_neutralize`): 14 residues, e.g. the
  thioesters SM099/SM117/SM221 (`[S+]=C-[O-]` -> `S-C=O`), 2MSU/BGAU/BGCU, MTDO, SM054, SM152, SM212.
- Still unusual after M3 (left as is): B1MA (`=[N+]=` cumulene in 1-methyladenine), PEGM (polymer fragment with open
  valences: `[CH2-][O+]=CH2`), SM212 (reference PG2 although the neutral form has one O-).
- 23 residues are disconnected (ions, water, multi-fragment) and 4 fail bond-order assignment; all are skipped.

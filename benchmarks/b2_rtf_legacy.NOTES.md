# B2: legacy RTF model compounds (no SMILES in the RTF)

Built by `src/open_cgenff/benchmarks/b2.py` from `top_all36_cgenff.rtf`; only residue names are stored here
(`b2_rtf_legacy.split.json`, 728 train / 182 test, seed 42). B2-test is sealed: aggregate metrics only.

Input quality (bond orders come from RDKit DetermineBondOrders on connectivity + net charge, not from the RTF);
issues seen on B2-train while checking perception in M2 (handle in the builder, never by editing reference types):
- PNTM (pentamidine-like bis-amidinium): one ring gets a quinoid structure with an `[O+]=` oxonium ether.
- SM033 (ethyl azide): azide drawn as `N=[N+]=[N-]`, pilot types differ by perception.
- 23 residues are disconnected (ions, water, multi-fragment) and 4 fail bond-order assignment; all are skipped.

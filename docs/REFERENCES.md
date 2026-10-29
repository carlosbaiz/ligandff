# References

PDFs are copyrighted and are not part of the repository. Carlos's copies are in the parent folder
(`CGenFF-compatible/`); file names below.

| id | citation | file | provides |
|---|---|---|---|
| P-I | Vanommeslaeghe K, MacKerell AD Jr. Automation of the CHARMM General Force Field (CGenFF) I: bond perception and atom typing. *J Chem Inf Model* 2012, 52, 3144–3154. doi:10.1021/ci300363c | `ci300363c.pdf`; SI `ci300363c_si_001.zip` (126 special-case test molecules, mol2, no reference types) | ring perception, resonance resolution + score, aromaticity rules, rule language, Fig. 5 rule excerpts |
| P-II | Vanommeslaeghe K, Raman EP, MacKerell AD Jr. Automation of the CHARMM General Force Field (CGenFF) II: assignment of bonded parameters and partial atomic charges. *J Chem Inf Model* 2012, 52, 3155–3168. doi:10.1021/ci3003649 | `ci3003649.pdf`; SI `ci3003649_si_001.pdf` (Table S1 increment constraints; Tables S2–S4 example .str outputs, program 0.9.x) | TPS, penalty-matrix rule format, bonded/nonbonded matrices, bond groups, charge-increment scheme, per-charge penalty eq 4, increment training |
| P-v5 | Croitoru A, Kumar A, Lambry J-C, Lee J, Sharif S, Yu W, MacKerell AD Jr, Aleksandrov A. Increasing the accuracy and robustness of the CHARMM General Force Field with an expanded training set. *J Chem Theory Comput* 2025, 21, 3044–3065. doi:10.1021/acs.jctc.5c00046 | `ct5c00046.pdf`; SI `ct5c00046_si_001.pdf` (Table S10: program 2.5.1 and 5.0 outputs for 4 molecules) | CGenFF 5.0: 1,390 new molecules, retrained increments, fused 4-ring rule, retraining details |
| P-X | Soteras Gutiérrez I, Lin F-Y, Vanommeslaeghe K, Lemkul JA, Armacost KA, Brooks CL III, MacKerell AD Jr. Parametrization of halogen bonds in the CHARMM general force field: improved treatment of ligand–protein interactions. *Bioorg Med Chem* 2016, 24, 4812–4825. doi:10.1016/j.bmc.2016.06.034 | `1-s2.0-S0968089616304576-main.pdf`, SI `…-mmc1.pdf` | halogen lone pairs: charge +0.05 e, distances, LJ |
| P-Orr | Orr AA, Sharif S, Wang J, MacKerell AD Jr. Preserving the integrity of empirical force fields. *J Chem Inf Model* 2022, 62, 3825–3831. doi:10.1021/acs.jcim.2c00615 | `ci2c00615.pdf`, SI `ci2c00615_si_001.pdf` | fidelity failures of MATCH; naming request; per-ligand type/charge tables |
| P-K | Knight JL, Yesselman JD, Brooks CL III. Assessing the quality of absolute hydration free energies among CHARMM-compatible ligand parameterization schemes. *J Comput Chem* 2013, 34, 893–903. doi:10.1002/jcc.23199 | `J Comput Chem - 2013 - Knight - ….pdf` | MATCH learned rules (1–2 bond environments) and their limits |
| CGenFF | Vanommeslaeghe K, et al. CHARMM General Force Field (CGenFF). *J Comput Chem* 2010, 31, 671–690. doi:10.1002/jcc.21367 | not provided | force-field design, parametrization protocol |

Data: https://github.com/Lemkul-Lab/toppar_c36_releases (CHARMM36 + CGenFF 5.0 files) and
https://github.com/mackerell-lab/cgenff5 (CGenFF 5.0 paper data); commits pinned in `scripts/fetch_data.sh`.

# Verification and experiment programs

All commands below are intended to be run from the repository root.

| Program | Purpose |
|---|---|
| `run_verifications.py` | Orchestrates the reproducibility checks and optional experiments. |
| `verify_gamma_criterion.py` | Exact/rational verification of the criterion with `alpha_c = 0.1547`, including improved first-cycle controls and tail bounds. |
| `compare_cordero2015.py` | Quantitative comparison with the semilocal criterion of Cordero et al. (2015). |
| `verify_exact_obstruction.py` | Verifies the counterexample to universal exact residual majorization. |
| `verify_algebra.py` | Checks algebraic identities used in the residual decomposition. |
| `symbolic_majorant_verification.py` | Symbolic checks for the scalar majorant and related expansions. |
| `verify_audit.py` | High-precision audit calculations used in the review of previous semilocal analyses. |
| `independent_verification_tjm2022.py` | Independent checks related to the 2022 analysis discussed in S4. |
| `check_additional_2022.py` | Additional exact checks for the 2022 material. |
| `reproduce_experiments.py` | Reproduces the five `N=30` experiments and the Chandrasekhar/Bratu refinement study. |
| `gamma_controls.py` | Shared exact-control definitions used by the criterion checks. |

Run `python code/run_verifications.py --help` for the main entry point.

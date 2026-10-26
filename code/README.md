# Programs

Run commands from the repository root. Exact checks and numerical reproductions are separate types of evidence. The complete derivations are in the supplement.

| Program | Section | Arithmetic and role |
|---|---|---|
| `run_verifications.py` | S1–S5 | Entry point; checks exit statuses and writes a classified execution summary. |
| `gamma_controls.py` | S3 | Shared rational controls, domain checks, outward rounding and tail functions. |
| `verify_gamma_criterion.py` | S3 | Exact first cycle, rational states and uniform tail inequalities for `1547/10000`. |
| `compare_cordero2015.py` | S3 | Exact comparisons for the same quadratic operator, norm and initial point. |
| `verify_exact_obstruction.py` | S1 | Fraction-based counterexample to exact scalar residual majorization. |
| `verify_lotfi_orders.py` | S1 | Exact power-series coefficients for M6/M9/M12; frozen matrices and ordered products. |
| `check_lotfi_numerically.py` | S1 | Independent 200-digit scaled-step checks; not a proof or interval validation. |
| `verify_algebra.py` | S1, S3 | Exact SymPy remainder, residual and tail-polynomial identities. |
| `symbolic_majorant_verification.py` | S2 | Exact rational factorizations and asymptotic constants; interval signs are proved analytically. |
| `verify_published_results.py` | S1, S4 | Exact counterexamples and separately labelled numerical reproductions for Cordero and Singh. |
| `verify_maroju2022.py` | S4 | Symbolic checks of selected estimates and the coupled local expansion. |
| `verify_maroju_bounds.py` | S4 | Exact power bounds and weak-regularity example; numerical radii labelled separately. |
| `reproduce_experiments.py` | S5 | Eleven Nyström runs, complete iterates and 401-point profiles. |
| `verify_experiment_data.py` | S5 | Source hashes, record completeness, increments, ACOC and printed table values. |
| `check_integrity.py` | Package | Byte-level checksums, independent of mathematical and numerical validation. |

Use `python code/run_verifications.py --help` for the orchestration options. `--rational-only` requires no third-party package. The full suite requires the pinned dependencies in `requirements.txt`.

The algorithms compute exact method stages in the mathematical model. Comparison defects in the proofs are not user-prescribed tolerances or modifications of the iteration.

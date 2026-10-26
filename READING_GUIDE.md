# Reading guide

This guide identifies the supporting material for **Semilocal convergence beyond exact majorization for an efficient three-step method under gamma conditions**. References in the first column are to the article. References beginning with `S` are to the standalone supplement; its equation and bibliography numbers are independent.

| Article | Supporting material | Program or data |
|---|---|---|
| Proposition 2.3; equations (6)–(8) | S1.1–S1.3: local expansion, quartic term and directional example | `code/verify_published_results.py`; `code/verify_algebra.py` |
| Proposition 2.4; Section 2.2 | S1.3 and the invariant-line reduction in S5.1, Problem 2 | Polynomial Hammerstein records in `results/experiments/` |
| Proposition 2.6 and Corollary 2.7 | S1.4: Q-order/R-order distinction and the lower-order obstruction | Analytical argument in the supplement |
| Example 2.9 | S1.4: failure of exact residual majorization under valid gamma conditions | `code/verify_exact_obstruction.py`; `results/exact_obstruction.json` |
| Introduction: the M6, M9 and M12 constructions of Lotfi et al. | S1.5, Proposition S1.1: coupled-field expansions and the scope of the counterexample | `code/verify_lotfi_orders.py`; `results/lotfi_orders.json` |
| Proposition 3.1 and Theorem 3.2; equations (13)–(18) | S2: scalar maps, polynomial signs and asymptotic constants | `code/symbolic_majorant_verification.py` |
| Corollary 3.3 | S2.6: summability of scalar controls | Analytical argument in the supplement |
| Remark 4.1 and Lemmas 4.3–4.5 | S3.1: gamma/Lipschitz condition and integrated remainder bounds | `code/verify_algebra.py` |
| Lemmas 4.6–4.9; equations (27)–(32) | S3.2: residual identities, recurrence and domain checks | `code/gamma_controls.py`; `code/verify_algebra.py` |
| Lemmas 4.10–4.12, Theorem 4.13 and Corollary 4.14 | S3.3 and S6.4: bounded controls, summable defects and error bounds | Exact criterion checks and analytical proofs |
| Lemma 4.15; equations (40)–(45) | S3.4: finite verification and the uniform tail induction | `code/verify_gamma_criterion.py` |
| Lemma 5.1, Corollary 5.2 and Table 1; equations (46)–(50) | S3.3–S3.5: first cycle, rational states, tail and localization | `results/criterion_alpha1547.json`, `.txt` |
| Section 5.2 and Table 2; equations (51)–(53) | S3.6: same operator, norm and initial point | `code/compare_cordero2015.py`; `results/comparison_2015.json` |
| Introduction and Remark 5.3: published semilocal analyses | S4: distinct scopes of false estimates, repairable arguments and refuted conclusions | `code/verify_published_results.py`, `code/verify_maroju2022.py`, `code/verify_maroju_bounds.py` |
| Section 6.2, Problems 1–5 | S5.1: operators, derivatives, domains and continuous gamma hypotheses | Parameters and models in `code/reproduce_experiments.py` |
| Table 3 | S5.4: all six states and ACOC histories for five problems at N=30 | `results/experiments/five_problems_N30.csv` and individual experiment folders |
| Table 4 | S5.3 and S5.5: reconstruction and refinement | `results/experiments/refinement.csv`, `profiles_grid401.csv`, `Bratu_reference.json` |
| Introduction: additive majorizations | S6: Argyros and Ferreira–Svaiter comparisons | Analytical source comparisons, with separate notation for each framework |

## Distinctions needed when reading the data

A nonzero one-step homogeneous term rules out a uniform local bound of higher degree. An exact Q-order for a particular orbit requires additional trajectory-level information.

The field `F(a,b)=(a+b**2,b+a**2)` demonstrates the failure of exact gamma residual majorization for the three-step method. The coupled field `F(a,b)=(a+a*b,b+a**2+a*b)` is used for the simultaneous Lotfi order counterexamples. On the path `(t,2*t)`, the first field does not refute sixth order for M6.

The hatted controls use polynomial bounds only in the initial cycle. They do not define a different algorithm. Their global admissibility is checked before applying the semilocal conclusion.

S4 concerns published mathematical statements, not numerical failure of every example. In particular, a failed sufficient condition is not evidence of divergence.

The five continuous problems live in `C[0,1]`. The computed nodal vectors solve Nyström discretizations. The analysis of the continuous operators and the experimental observations are different claims.

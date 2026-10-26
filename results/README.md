# Results

## Exact rational and symbolic material

- `criterion_alpha1547.json` and `.txt`: first-cycle values, rational states and tail inequalities.
- `comparison_2015.json`: applicability intervals and location bounds for the fixed quadratic family.
- `exact_obstruction.json`: exact-majorization counterexample.
- `lotfi_orders.json`: all retained rational Taylor coefficients and leading terms for two quadratic fields.
- `algebra_verification.json`: exact symbolic identities.

## Published-results checks

`literature_checks/published_results.json` separates exact assertions from multiprecision reproductions. `literature_checks/cordero2015_iterations.csv` contains the six-cycle eight-node reproduction. `maroju_bounds.json` distinguishes exact inequalities from numerically computed radii. These are mathematical checks of cited publications, with their scope explained in S4.

`lotfi_numerical.json` contains direct high-precision evaluations of the Lotfi coupled example. They are not used to prove the exact Taylor orders.

## Numerical experiments

`experiments/` contains eleven experiment folders and the aggregate files:

- `five_problems_N30.csv`: six-state histories for the five problems at N=30.
- `refinement.csv`: eight Chandrasekhar/Bratu profile comparisons.
- `profiles_grid401.csv`: reconstructed profiles on the common grid.
- `Bratu_reference.json`: the lower-amplitude reference parameter and equation.
- `RUN_STATUS.json`: protocol, source-digest, increment, ACOC and rounding checks.

All reported numerical residuals and profile differences are multiprecision floating-point quantities, not interval certificates or quadrature-error bounds.

## Execution records

`verification/summary.json` identifies the arithmetic and exit status of the executed verification commands. Its accompanying logs give their output. `experiment_execution.json` records the separate eleven-experiment execution. `supplement_build.json` records the supplement's compilation checks. These records describe the distributed computation, not a performance benchmark.

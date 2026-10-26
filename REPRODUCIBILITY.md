# Reproducibility guide

## Environment

The supplied executions use Python 3.13.5 on Linux x86_64, `mpmath==1.3.0` (Python backend) and `sympy==1.14.0`. Python 3.10 or later is supported by the source syntax. Install the pinned dependencies with:

```bash
python -m pip install -r requirements.txt
```

The exact rational programs use the standard-library `fractions.Fraction`. No floating-point comparison decides their acceptance inequalities. SymPy programs check identities symbolically; the polynomial sign proofs on intervals remain analytical arguments in S2. Multiprecision results are explicitly distinguished from these exact checks.

Run commands from the repository root. Use ordinary Python execution, without `-O` or `PYTHONOPTIMIZE`, so the assertion checks remain active.

## Check the distributed files

```bash
python code/check_integrity.py
```

This checks every entry of `SHA256SUMS.txt`. Check it before regeneration: execution times, platform strings and PDF build metadata can change when outputs are recomputed. A different checksum for a regenerated file is not itself a mathematical discrepancy. The numerical data have separate protocol, coefficient, increment and rounding checks.

## Rational verification

```bash
python code/run_verifications.py --rational-only
```

This executes:

1. The first-cycle bounds, two rational cycles and the uniform tail inequalities for `alpha_c=1547/10000`.
2. The comparison with Cordero et al. (2015), for the same operator, norm and initial point.
3. The exact-majorization obstruction at `h=1/1000`.
4. The six one-step Lotfi expansions for the two quadratic fields, using exact truncated power series.

The main certificates are `results/criterion_alpha1547.json`, `results/comparison_2015.json`, `results/exact_obstruction.json` and `results/lotfi_orders.json`. Decimal presentations of rational values are not used to determine signs.

## Full non-experimental verification

```bash
python code/run_verifications.py
python -m unittest discover -s tests -v
```

This also checks integrated remainders, the tail polynomial, scalar factorizations, published estimates, counterexamples and numerical reproductions described in S4. The Cordero eight-node calculation uses 2000 digits. Its output is labelled as numerical rather than interval arithmetic. Direct 200-digit evaluations independently check the Lotfi scaled coefficients; the proof is the exact expansion, not these samples.

Each subprocess writes a log in `results/verification/`. The runner stops on the first nonzero exit status. `summary.json` records the command, arithmetic category, source digest, elapsed time and exit status of each requested check. `--rational-only` replaces this summary with the rational subset; running the default command restores a full non-experimental summary.

## Eleven Nyström experiments

```bash
python code/run_verifications.py --experiments
```

This runs the full verification suite, all five problems at N=30, and the Chandrasekhar and Bratu problems at N=20, 40 and 60. The experiment program can also be used directly:

```bash
python code/reproduce_experiments.py --problem all --N 30 --out results/experiments
python code/reproduce_experiments.py --problem chandrasekhar --N 20 40 60 --out results/experiments
python code/reproduce_experiments.py --problem bratu --N 20 40 60 --out results/experiments
python code/reproduce_experiments.py --summarize --out results/experiments
python code/verify_experiment_data.py
```

Use a separate output directory for exploratory node counts, precisions or numbers of cycles; do not mix such runs with the eleven-record dataset used by `--summarize` and the table checks.

The protocol is fixed: Gauss–Legendre quadrature on `[0,1]`, 1000 decimal digits for nodes and iteration, five cycles, maximum norm of nodal vectors and a single pivoted LU factorization reused for three solves per cycle. No lower-precision quadrature nodes or binary double-precision parameters are introduced. Profile reconstruction uses 100 digits on the 401-point grid `{i/400: 0<=i<=400}`.

Each experiment folder contains:

| File | Contents |
|---|---|
| `result.json` | Parameters, six-state history, nodes, weights, solution, numerical moment/linear-solve residuals, software/platform and source digest. |
| `complete_iterates.json` | Full-precision main iterates and all intermediate stages. |
| `history.csv` | Residual, increment, ACOC and precision-status flags at each state. |
| `nodes_weights_solution.csv` | Nodes, quadrature weights and final nodal vector. |

Aggregate files are `five_problems_N30.csv`, `refinement.csv`, `profiles_grid401.csv` and `Bratu_reference.json`. The data checker verifies the eleven records, source digests, 55 increments, 33 ACOC values, the printed Table 3 and Table 4 roundings, and the 3208 profile rows. It writes `results/experiments/RUN_STATUS.json`.

## Precision and interpretation

The ACOC uses three consecutive nonzero increments and is omitted when its denominator vanishes or an increment reaches the stated precision floor. A residual near the precision floor is flagged, not converted into an exact bound of the form `10**(-dps)`.

The Bratu reference is the lower-amplitude positive solution

$$u(s)=2\log\frac{\cosh(a/2)}{\cosh(a(s-1/2))},\qquad
\frac{2a^2}{\cosh^2(a/2)}=2,\quad 0<a<2.$$

For Chandrasekhar, the N=60 reconstruction is the numerical reference. Grid differences are not continuous sup-norm bounds. Discrete residuals are not quadrature-error bounds, and none of the multiprecision results are interval computations. Execution times are recorded for provenance, not as a performance comparison.

The polynomial Hammerstein equation is the model identified in Cárdenas–Castro–Sierra, Example 2, equation (3.3), p. 4617. The reciprocal Hammerstein equation uses the kernel and reciprocal nonlinearity of Amat–Busquier–Gutiérrez, equation (6), p. 2941, with coefficient `mu=0.1` and the protocol specified here.

## Build the supplement

A TeX distribution needs `pdflatex`, the standard AMS packages, `mathtools`, `lmodern`, `microtype`, `booktabs`, `array`, `longtable`, `geometry`, `xcolor`, `xurl`, `hyperref`, `fancyhdr`, `etoolbox`, `babel` and the corresponding fonts.

```bash
python build_supplement.py
```

The builder uses three passes in a temporary directory, checks the final log for unresolved references and overfull boxes, and copies only the PDF into the repository root. The source has no dependency on the article. No article source, external labels, proprietary fonts or publisher class files are required.

The programs and mathematical proofs are complementary. Finite checks do not replace the uniform tail induction or the analytical interval-sign proofs.

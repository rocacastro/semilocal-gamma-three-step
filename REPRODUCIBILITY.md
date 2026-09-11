# Reproducibility guide

## 1. Software environment

The supplied runs were produced with Python 3.13.5 on Linux x86_64 using:

- `mpmath==1.3.0`
- `sympy==1.14.0`

The programs are written for Python 3.10 or later. No internet connection is required after the dependencies are installed.

## 2. Exact and symbolic checks

From the repository root:

```bash
python code/run_verifications.py --rational-only
```

This checks the explicit gamma criterion with

\[
\alpha_c=1547/10000=0.1547,
\]

including the improved first-cycle controls, and the rational comparison with the 2015 criterion.

The broader check

```bash
python code/run_verifications.py
```

also runs the exact obstruction example, the algebraic identities, the symbolic majorant checks, and the audits documented in Supplementary Sections S1–S4 and S6.

The combined logs are written to `results/combined_run/`.

## 3. Numerical experiments in S5

To reproduce all eleven experiments:

```bash
python code/run_verifications.py --experiments
```

This runs:

- all five problems with `N=30`;
- Chandrasekhar with `N=20,40,60` in addition to `N=30`;
- Bratu with `N=20,40,60` in addition to `N=30`;
- the final summary and 401-point profile comparison.

The numerical protocol uses Gauss–Legendre Nyström discretization on `[0,1]`, 1000 decimal digits for iterations, five cycles of the three-step method, one LU factorization per cycle with three solves, the maximum norm for nodal vectors, and 100 digits for reconstructed profiles on the 401-point grid.

Main summary files:

- `results/experiments/five_problems_N30.csv`
- `results/experiments/refinement.csv`
- `results/experiments/profiles_grid401.csv`
- `results/experiments/Bratu_reference.json`

Each individual experiment directory contains `result.json`, `history.csv`, `complete_iterates.json`, and `nodes_weights_solution.csv`.

## 4. Interpretation

The continuous semilocal constants are established analytically. The Nyström runs are numerical illustrations of the corresponding discretized problems. The two levels should not be conflated.

The approximate computational order of convergence is computed from successive increments. Values affected by zero increments, a zero denominator, or proximity to the working-precision floor are not used.

The Bratu reference profile uses the small exact branch

\[
u(s)=2\log\frac{\cosh(a/2)}{\cosh(a(s-1/2))},
\qquad
\frac{2a^2}{\cosh^2(a/2)}=2.
\]

## 5. Building the Supplementary Material

Run:

```bash
python build_supplement.py
```

or compile manually:

```bash
pdflatex -interaction=nonstopmode -halt-on-error supplementary_material_gamma.tex
pdflatex -interaction=nonstopmode -halt-on-error supplementary_material_gamma.tex
pdflatex -interaction=nonstopmode -halt-on-error supplementary_material_gamma.tex
```

The Supplementary Material is self-contained. Generated `.aux`, `.log`, `.out`, and related TeX files are excluded from version control.

## 6. Release procedure for Zenodo

Recommended sequence:

1. Upload this repository to GitHub.
2. Confirm that the repository files and results are unchanged.
3. Create a GitHub release tagged `v1.0.0`.
4. Archive that release in Zenodo.
5. Add the assigned DOI to `CITATION.cff` and `.zenodo.json` if desired for a subsequent repository release.

The archived Zenodo release itself provides the persistent identifier for the exact repository version.

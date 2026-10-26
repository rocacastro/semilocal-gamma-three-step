# Semilocal convergence beyond exact majorization

Supplementary material and reproducibility package for

**Semilocal convergence beyond exact majorization for an efficient three-step method under gamma conditions**

Rodrigo Castro Marín · Lía Paz Guzmán Pinto · Luis Fernando Echeverri Delgado · Óscar Antonio Restrepo Gutiérrez · Gerardo Andrés Honorato Gutiérrez

This research compendium provides the mathematical derivations, verification programs and numerical data supporting the article. Start with [`supplementary_material_gamma.pdf`](supplementary_material_gamma.pdf) and the [reading guide](READING_GUIDE.md), which links the article's results to the supplementary sections and executable checks.

The sufficient condition established by the exact rational checks is

$$0<\beta\gamma\leq\frac{1547}{10000}=0.1547.$$

## Mathematical contents

| Section | Contents |
|---|---|
| **S1** | Local error expansion with ordered multilinear compositions; the obstruction to exact majorization; an exact gamma counterexample; one-step expansions for the M6, M9 and M12 maps of Lotfi et al. |
| **S2** | Complete factorizations, sign arguments, interlacing, asymptotic constants and summability for the scalar reference sequences. |
| **S3** | Effective-radius controls, polynomial first-cycle bounds, exact rational verification of the criterion, uniform tail control and comparison with Cordero et al. (2015). |
| **S4** | Mathematical verification of specific results and estimates in the semilocal analyses of Cordero et al. (2015), Singh et al. (2016) and Maroju et al. (2022). |
| **S5** | Continuous hypotheses, derivatives, eleven Nyström experiments, complete iterates and profile refinements for five integral equations. |
| **S6** | Comparison with Argyros and Ferreira–Svaiter; the meaning and logical role of summable comparison defects. |

The Lotfi checks use the coupled quadratic field `F(a,b)=(a+a*b, b+a**2+a*b)` and retain the frozen matrices and the multiplication order of the published algorithms. Nonzero coefficients of degrees 5, 7 and 9 exclude the corresponding uniform local bounds of degrees 6, 9 and 12. They do **not** assign an exact Q-order to every trajectory.

## Quick start

Python 3.10 or later is required. The supplied execution uses Python 3.13.5, `mpmath==1.3.0` and `sympy==1.14.0`.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell:
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Run the rational checks, including the criterion, the comparison, the exact-majorization counterexample and the Lotfi expansions:

```bash
python code/run_verifications.py --rational-only
```

Run the full non-experimental suite:

```bash
python code/run_verifications.py
python -m unittest discover -s tests -v
```

Recompute the eleven numerical experiments and check their records:

```bash
python code/run_verifications.py --experiments
```

The experiments use 1000 decimal digits, five cycles and one LU factorization per cycle. Numerical verification of the stored data can also be run separately:

```bash
python code/verify_experiment_data.py
```

Do not invoke the programs with Python's `-O` option: assertions are part of the mathematical checks. No network access is required after installation of the dependencies.

## Files and outputs

| Location | Purpose |
|---|---|
| `supplementary_material_gamma.tex`, `.pdf` | Standalone mathematical supplement, S1–S6. |
| `code/` | Rational, symbolic, multiprecision and integrity-checking programs. |
| `tests/` | Regression tests and rounded numerical values printed in the article. |
| `results/` | Exact certificates, labelled numerical reproductions, complete iterates and execution records. |
| `READING_GUIDE.md` | Correspondence between article results, supplement and programs. |
| `REPRODUCIBILITY.md` | Commands, output interpretation, precision policy and data layout. |
| `SHA256SUMS.txt` | SHA-256 digests of the distributed files. |
| `CITATION.cff`, `.zenodo.json`, `VERSION` | Citation and package metadata. |

The [program index](code/README.md) and [results index](results/README.md) describe individual files. To check the distributed bytes before rerunning programs:

```bash
python code/check_integrity.py
```

A fresh execution may change elapsed times and environment metadata; such regenerated files need not have the same byte-level digest as the distributed records. Mathematical checks and comparisons with printed values are separate from checksum verification.

## Compile the supplement

Install a TeX distribution with `pdflatex` and the packages listed in the [reproducibility guide](REPRODUCIBILITY.md), then run:

```bash
python build_supplement.py
```

The source compiles without the article's source, class, labels or auxiliary files. Mathematical references in the supplement have their own `S`-numbering; its bibliography is numbered independently.

## Scope of the evidence

The rational checks use `fractions.Fraction`; the symbolic checks use exact SymPy expressions. The multiprecision calculations are **not** interval computations. Residuals of the discretized equations do not bound quadrature error, and maxima on a 401-point grid are not certified continuous sup-norm errors. The effective-radius bound is no larger than the bound expressed through comparison defects; the use of a positive part alone is not a quantitative improvement.

The Python checks support the analytical proofs in S1–S6; finite execution is not a replacement for the uniform tail argument.

## Citation and licensing

Use the author and title information in [`CITATION.cff`](CITATION.cff). The repository address is `https://github.com/rocacastro/semilocal-gamma-three-step`.

File-level licensing is described in [`LICENSE.md`](LICENSE.md). The MIT software license does not automatically license the scholarly text or the numerical data.

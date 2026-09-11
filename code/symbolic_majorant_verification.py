"""Symbolic verifications for Supplementary Material S2.

The program uses SymPy to check the two main rational factorizations
and the evaluations yielding the asymptotic constant of T.
It does not replace the sign proofs on the admissible interval.
"""
from __future__ import annotations

import sympy as sp

x, r = sp.symbols("x r")
alpha = r * (1 - 2 * r) / (1 - r)

G = sp.cancel(alpha - x + x**2 / (1 - x))
Gp = sp.diff(G, x)
S = sp.cancel(x - G / Gp)
GS = sp.cancel((alpha - S + S**2 / (1 - S)))
V = sp.cancel(S - 5 * GS / Gp)
GV = sp.cancel(alpha - V + V**2 / (1 - V))
T = sp.cancel(S - (9 * GS + GV) / (5 * Gp))

D = 2*r**2*x - 2*r**2 - 4*r*x + 2*r + 3*x - 1
R = (
    32*r**4*x**4 - 128*r**4*x**3 + 216*r**4*x**2 - 176*r**4*x + 56*r**4
    - 128*r**3*x**4 + 428*r**3*x**3 - 612*r**3*x**2 + 414*r**3*x - 102*r**3
    + 204*r**2*x**4 - 564*r**2*x**3 + 708*r**2*x**2 - 426*r**2*x + 87*r**2
    - 152*r*x**4 + 314*r*x**3 - 318*r*x**2 + 169*r*x - 31*r
    + 52*x**4 - 82*x**3 + 60*x**2 - 25*x + 4
)
Q = (
    8*r**4*x**4 - 32*r**4*x**3 + 24*r**4*x**2 + 16*r**4*x - 16*r**4
    - 32*r**3*x**4 + 152*r**3*x**3 - 168*r**3*x**2 + 36*r**3*x + 12*r**3
    + 36*r**2*x**4 - 216*r**2*x**3 + 252*r**2*x**2 - 84*r**2*x + 3*r**2
    - 8*r*x**4 + 116*r*x**3 - 132*r*x**2 + 46*r*x - 4*r
    - 2*x**4 - 28*x**3 + 30*x**2 - 10*x + 1
)
P = (
    8*r**3*x**5 - 40*r**3*x**4 + 56*r**3*x**3 - 8*r**3*x**2 - 32*r**3*x + 16*r**3
    + 24*r**2*x**6 - 168*r**2*x**5 + 424*r**2*x**4 - 424*r**2*x**3
    + 100*r**2*x**2 + 72*r**2*x - 28*r**2
    - 48*r*x**6 + 276*r*x**5 - 548*r*x**4 + 382*r*x**3 - 2*r*x**2 - 79*r*x + 19*r
    + 20*x**6 - 92*x**5 + 122*x**4 - 6*x**3 - 63*x**2 + 33*x - 5
)

Vprime_rhs = 2*(r-x)*(1-2*r-2*x+2*r*x)*R / (
    (r-1)*(2*x**2-4*x+1)**3*D**2
)
Tminus_r_rhs = 2*(r-x)**5*P / (
    (r-1)*(2*x**2-4*x+1)**3*D*Q
)


def rational_identity(lhs: sp.Expr, rhs: sp.Expr, name: str) -> None:
    num = sp.together(lhs-rhs).as_numer_denom()[0]
    ok = sp.Poly(sp.expand(num), x, r).is_zero
    print(f"{name}: {ok}")
    if not ok:
        raise AssertionError(f"Failure in {name}")


rational_identity(sp.diff(V, x), Vprime_rhs, "Factorization of V'(x)")
rational_identity(T-r, Tminus_r_rhs, "Factorization of T(x)-r")

d = 2*r**2 - 4*r + 1
checks = {
    "D_r(r)": sp.factor(D.subs(x, r) - (r-1)*d),
    "Q_r(r)": sp.factor(Q.subs(x, r) - (1-r)**2*d**3),
    "P_r(r)": sp.factor(P.subs(x, r) - (2*r-5)*(2*r+1)*d**3),
}
for name, value in checks.items():
    ok = value == 0
    print(f"{name}: {ok}")
    if not ok:
        raise AssertionError(f"Failure in {name}: {value}")

# The constant is obtained directly from the factorization of T-r.
constant_from_factor = sp.factor(
    -2 * P.subs(x, r)
    / ((r-1) * d**3 * D.subs(x, r) * Q.subs(x, r))
)
constant_expected = 2*(1+2*r)*(5-2*r) / ((1-r)**4*d**4)
constant_ok = sp.factor(constant_from_factor-constant_expected) == 0
print(f"Normalized asymptotic constant of T: {constant_ok}")
if not constant_ok:
    raise AssertionError(sp.factor(constant_from_factor-constant_expected))

print("All symbolic identities were verified.")

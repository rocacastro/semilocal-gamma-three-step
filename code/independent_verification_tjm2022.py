#!/usr/bin/env python3
"""Exact checks of two counterexamples for scheme (2.1).

Source examined: P. Maroju, R. Behl and S. S. Motsa,
Thai Journal of Mathematics 20(1) (2022), 21-33.

Requirement: SymPy (python -m pip install sympy).
Run: python independent_verification_tjm2022.py

Logical scope:
- The vector counterexample excludes a LOCAL UNIFORM estimate O(||e||^5).
  It does not assert that all trajectories have the same Q-order or R-order.
- The first scalar counterexample refutes estimate (2.6), even when phi is
  extended to [0,infinity) compatibly with the printed C6.
- The check of (2.7) uses phi(s)=s: C6 and its domain must be repaired.
  This is NOT a counterexample literally satisfying all the incompatible
  published hypotheses, nor does it prove divergence of the method.
"""
from __future__ import annotations

try:
    import sympy as sp
except ImportError as exc:
    raise SystemExit("Install the dependency with: python -m pip install sympy") from exc


def verify() -> None:
    t = sp.symbols("t", real=True)
    a = sp.symbols("a", positive=True)
    degree = 6

    def truncate(vector: sp.Matrix) -> sp.Matrix:
        return vector.applyfunc(
            lambda c: sp.series(c, t, 0, degree).removeO().expand()
        )

    def field(vector: sp.Matrix) -> sp.Matrix:
        u, v = vector
        return sp.Matrix([u + v**2, v + u**2])

    # F(u,v)=(u+v^2,v+u^2), x=(t,2t), F'(0)=I.
    x = sp.Matrix([t, 2*t])
    jacobian = sp.Matrix([[1, 4*t], [2*t, 1]])
    inverse = jacobian.inv()
    y = truncate(x - inverse*field(x))
    z = truncate(y - 5*inverse*truncate(field(y)))
    published_step = truncate(z - inverse*(-16*truncate(field(y)) + truncate(field(z)))/5)
    equivalent_step = truncate(y - inverse*(9*truncate(field(y)) + truncate(field(z)))/5)
    expected = sp.Matrix([28*t**4 + 224*t**5, -56*t**4 + 70*t**5])
    assert published_step == equivalent_step == expected
    assert sp.factor(jacobian.det()) == 1 - 8*t**2
    print("1. Equivalence of the steps and expansion through degree five: OK.")
    print("   T(t,2t) =", expected.T, "+ O(t^6)")
    print("   Limit ||T(t,2t)||_inf / ||(t,2t)||_inf^4 =", sp.Rational(56,16))

    # Independent identification of the fourth-degree homogeneous term.
    u, v = sp.symbols("u v", real=True)
    def bilinear(p: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
        return sp.Matrix([p[1]*q[1], p[0]*q[0]])
    e = sp.Matrix([u,v])
    b = bilinear(e,e)
    q4 = 4*(bilinear(e,bilinear(e,b))-bilinear(b,b))
    q4_expected = 4*(v**3-u**3)*sp.Matrix([u,-v])
    assert (q4-q4_expected).applyfunc(sp.expand) == sp.zeros(2,1)
    print("2. Homogeneous term Q4(e) = 4(v^3-u^3)(u,-v): OK.")

    # F_a(s)=s-1-a*s^2/2. At x0=0, Gamma0=beta=eta=1, a0=a.
    def fa(s):
        return s - 1 - a*s**2/2
    x0 = sp.S.Zero
    y0 = x0-fa(x0)
    z0 = sp.expand(y0-5*fa(y0))
    x1 = sp.expand(z0-(-16*fa(y0)+fa(z0))/5)
    assert y0 == 1 and z0 == 1+5*a/2
    assert x1 == 1+a/2+a**2/2+5*a**3/8

    # omega(t)=a*t and phi(s)=max(1,s): M=integral_0^1 phi=1.
    # This is an extension compatible with the published C6.
    lhs26 = a*z0**2
    rhs26 = a*(1+a)**2
    assert sp.factor(lhs26-rhs26) == 3*a**2*(7*a+4)/4
    a_test = sp.Rational(1,10)
    assert lhs26.subs(a,a_test) == sp.Rational(5,32)
    assert rhs26.subs(a,a_test) == sp.Rational(121,1000)
    assert lhs26.subs(a,a_test) > rhs26.subs(a,a_test)
    print("3. Bound (2.6), after extending the domain of phi: FALSE.")
    print("   For a=1/10, the left-hand side is 5/32 > 121/1000.")

    # Check of the supplied audit: phi(s)=s, M=1/2.
    # WARNING: this phi does NOT satisfy phi>=1 on [0,1].
    published_h = a/2+a**2/10+a**3/40
    assert sp.factor(x1-1-published_h) == a**2*(3*a+2)/5
    assert x1.subs(a,a_test) == sp.Rational(1689,1600)
    assert (1+published_h).subs(a,a_test) == sp.Rational(42041,40000)
    print("4. Bound (2.7) under the natural repair phi(s)=s: FALSE.")
    print("   x1-(1+h(a)) = a^2(3a+2)/5 > 0.")
    print("   For a=1/10: 1.055625 > 1.051025.")
    print("The exact checks completed successfully.")


if __name__ == "__main__":
    verify()

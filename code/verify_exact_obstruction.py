#!/usr/bin/env python3
"""Exact verification of the obstruction to residual majorization.

Uses only the standard library. Does not use floating-point numbers to
accept any inequality. The general asymptotic expansion is proved in S1;
here we verify a fixed datum satisfying the gamma conditions.

Run: python code/verify_exact_obstruction.py
"""
from __future__ import annotations
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
import json

Vector = tuple[Q, Q]
Matrix = tuple[tuple[Q, Q], tuple[Q, Q]]

def field(x: Vector) -> Vector:
    u, v = x
    return u + v*v, v + u*u

def inverse_jacobian(x: Vector) -> Matrix:
    u, v = x
    determinant = 1 - 4*u*v
    if determinant == 0:
        raise ArithmeticError('The derivative is singular.')
    return ((1/determinant, -2*v/determinant),
            (-2*u/determinant, 1/determinant))

def apply(a: Matrix, x: Vector) -> Vector:
    return (a[0][0]*x[0]+a[0][1]*x[1],
            a[1][0]*x[0]+a[1][1]*x[1])

def sub(x: Vector, y: Vector) -> Vector:
    return x[0]-y[0], x[1]-y[1]

def scale(a: Q, x: Vector) -> Vector:
    return a*x[0], a*x[1]

def norm(x: Vector) -> Q:
    return max(abs(x[0]), abs(x[1]))

def one_step(x: Vector) -> tuple[Vector, Vector, Vector]:
    a = inverse_jacobian(x)
    y = sub(x, apply(a, field(x)))
    z = sub(y, scale(Q(5), apply(a, field(y))))
    fy, fz = field(y), field(z)
    q = ((9*fy[0]+fz[0])/5, (9*fy[1]+fz[1])/5)
    return y, z, sub(y, apply(a, q))

def decimal(x: Q) -> str:
    with localcontext() as context:
        context.prec = 45
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def record(x: Q) -> dict[str, str]:
    return {'rational': str(x), 'decimal_for_display': decimal(x)}

def main() -> None:
    h = Q(1,1000)
    x0 = (h,2*h)
    beta = (2*h+h*h)/(1-4*h)
    gamma = 1/(1-4*h)
    alpha = beta*gamma
    a0 = inverse_jacobian(x0)
    matrix_norm = max(sum(abs(e) for e in row) for row in a0)
    y, z, x1 = one_step(x0)
    real_residual = norm(apply(a0, field(x1)))
    def g(t: Q) -> Q:
        return beta-t+gamma*t*t/(1-gamma*t)
    s0 = beta
    v0 = s0 + 5*g(s0)
    t1 = s0 + (9*g(s0)+g(v0))/5
    scalar_residual = g(t1)
    # Base hypotheses.  F'' is constant with norm 2; all differences vanish.
    assert 0 < h < Q(1,4)
    assert norm(field(x0)) == 2*h+h*h
    assert matrix_norm <= 1/(1-4*h)
    assert norm(apply(a0, field(x0))) <= beta
    assert 2*matrix_norm <= 2*gamma
    assert alpha == Q(667,330672)
    assert 0 < alpha < Q(1547,10000)
    assert (3-alpha)**2 > 8
    # First scalar cycle is inside the pole.
    assert 0 < s0 < 1/gamma and 0 < v0 < 1/gamma
    assert 0 < t1 < 1/gamma
    # Exact comparison; no rounded decimal participates in these assertions.
    assert real_residual > Q(5,10**11)
    assert 0 < scalar_residual < Q(4,10**13)
    assert real_residual > scalar_residual
    result = {
        'h':record(h), 'beta':record(beta), 'gamma':record(gamma),
        'alpha':record(alpha), 'x1':[record(x) for x in x1],
        'real_residual':record(real_residual),
        'scalar_residual':record(scalar_residual),
        'ratio':record(real_residual/scalar_residual),
        'status':'All comparisons were verified with exact fractions.'
    }
    folder = Path(__file__).resolve().parent.parent / 'results'
    folder.mkdir(exist_ok=True)
    (folder/'exact_obstruction.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(result['status'])
    print('Normalized residual: ',decimal(real_residual))
    print('Majorant residual: ',decimal(scalar_residual))
    print('Ratio: ',decimal(real_residual/scalar_residual))

if __name__ == '__main__':
    main()

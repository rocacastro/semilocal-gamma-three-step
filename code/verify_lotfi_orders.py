#!/usr/bin/env python3
"""Exact Taylor coefficients for the M6, M9 and M12 maps of Lotfi et al.

Source: DOI 10.1080/00207160.2014.946412, equations (5), (22), (24).
All operations are in Q[t]/(t**12), implemented with fractions.Fraction.
The frozen matrix H and the order of matrix products follow the source.
Nonzero one-step coefficients refute uniform local bounds; they do not
assign an exact Q-order to every complete iteration sequence.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

Q = Fraction
TERMS = 12
Series = list[Fraction]
Vector = list[Series]
Matrix = list[list[Series]]


def constant(value: int | Fraction = 0) -> Series:
    return [Q(value)] + [Q(0)] * (TERMS - 1)


def add(left: Series, right: Series) -> Series:
    return [a + b for a, b in zip(left, right)]


def negate(value: Series) -> Series:
    return [-coefficient for coefficient in value]


def subtract(left: Series, right: Series) -> Series:
    return add(left, negate(right))


def multiply(left: Series, right: Series) -> Series:
    result = constant()
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right[:TERMS - i]):
                if b:
                    result[i + j] += a * b
    return result


def reciprocal(value: Series) -> Series:
    """Invert a power series with nonzero constant coefficient."""
    if not value[0]:
        raise ZeroDivisionError('A power series inverse needs a nonzero constant term.')
    result = constant()
    result[0] = 1 / value[0]
    for degree in range(1, TERMS):
        result[degree] = -sum(
            (value[k] * result[degree - k] for k in range(1, degree + 1)), Q(0)
        ) / value[0]
    return result


def scale(value: Series, scalar: int | Fraction) -> Series:
    return [coefficient * scalar for coefficient in value]


def vector_subtract(left: Vector, right: Vector) -> Vector:
    return [subtract(a, b) for a, b in zip(left, right)]


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return [[add(left[i][j], right[i][j]) for j in range(2)] for i in range(2)]


def matrix_subtract(left: Matrix, right: Matrix) -> Matrix:
    return [[subtract(left[i][j], right[i][j]) for j in range(2)] for i in range(2)]


def matrix_scale(value: Matrix, scalar: int | Fraction) -> Matrix:
    return [[scale(entry, scalar) for entry in row] for row in value]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return [[add(multiply(left[i][0], right[0][j]),
                 multiply(left[i][1], right[1][j])) for j in range(2)] for i in range(2)]


def matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    return [add(multiply(matrix[i][0], vector[0]),
                multiply(matrix[i][1], vector[1])) for i in range(2)]


def determinant(matrix: Matrix) -> Series:
    return subtract(multiply(matrix[0][0], matrix[1][1]),
                    multiply(matrix[0][1], matrix[1][0]))


def matrix_inverse(matrix: Matrix) -> Matrix:
    inverse_determinant = reciprocal(determinant(matrix))
    adjugate = [[matrix[1][1], negate(matrix[0][1])],
                [negate(matrix[1][0]), matrix[0][0]]]
    return [[multiply(entry, inverse_determinant) for entry in row] for row in adjugate]


def identity() -> Matrix:
    return [[constant(1), constant()], [constant(), constant(1)]]


def leading_term(vector: Vector) -> tuple[int, Fraction, Fraction]:
    for degree in range(TERMS):
        if vector[0][degree] or vector[1][degree]:
            return degree, vector[0][degree], vector[1][degree]
    raise ArithmeticError('No nonzero coefficient occurs at the retained degrees.')


def analyze(coefficients: tuple[int, int, int, int, int, int]) -> dict:
    """Evaluate one cycle for a quadratic field with F(0)=0 and F'(0)=I.

    F(a,b)=(a+c0*a**2+c1*a*b+c2*b**2,
            b+c3*a**2+c4*a*b+c5*b**2).
    """
    a11, a12, a22, b11, b12, b22 = coefficients

    def field(vector: Vector) -> Vector:
        u, v = vector
        uu, uv, vv = multiply(u, u), multiply(u, v), multiply(v, v)
        return [add(u, add(scale(uu, a11), add(scale(uv, a12), scale(vv, a22)))),
                add(v, add(scale(uu, b11), add(scale(uv, b12), scale(vv, b22))))]

    def jacobian(vector: Vector) -> Matrix:
        u, v = vector
        return [
            [add(constant(1), add(scale(u, 2*a11), scale(v, a12))),
             add(scale(u, a12), scale(v, 2*a22))],
            [add(scale(u, 2*b11), scale(v, b12)),
             add(constant(1), add(scale(u, b12), scale(v, 2*b22)))]
        ]

    t = constant()
    t[1] = Q(1)
    x = [t, scale(t, 2)]
    jx = jacobian(x)
    inverse_jx = matrix_inverse(jx)
    y = vector_subtract(x, matrix_vector(inverse_jx, field(x)))
    jy = jacobian(y)
    sum_jacobians = matrix_add(jx, jy)
    z = vector_subtract(x, matrix_vector(matrix_scale(matrix_inverse(sum_jacobians), 2), field(x)))
    b = matrix_multiply(inverse_jx, jy)
    h = matrix_add(matrix_subtract(matrix_scale(identity(), Q(7, 2)), matrix_scale(b, 4)),
                   matrix_scale(matrix_multiply(b, b), Q(3, 2)))
    k = matrix_multiply(h, inverse_jx)
    outputs = {}
    point = z
    for name in ('M6', 'M9', 'M12'):
        point = vector_subtract(point, matrix_vector(k, field(point)))
        outputs[name] = point
    return {'outputs': outputs, 'z': z, 'I_minus_K': matrix_subtract(identity(), k),
            'jacobian_determinant': determinant(jx),
            'sum_jacobians_determinant': determinant(sum_jacobians)}


CASES = {
    'coupled_quadratic': {
        'field': '(a + a*b, b + a**2 + a*b)',
        'coefficients': (0, 1, 0, 1, 1, 0),
        'expected': [(5, Q(-3, 4), Q(3, 4)), (7, Q(9, 8), Q(-9, 8)),
                     (9, Q(-27, 16), Q(27, 16))],
    },
    'symmetric_quadratic': {
        'field': '(a + b**2, b + a**2)',
        'coefficients': (0, 0, 1, 1, 0, 0),
        'expected': [(6, Q(360), Q(-36)), (8, Q(-1512), Q(3024)),
                     (10, Q(12096), Q(-24192))],
    },
}


def verify() -> dict:
    report = {
        'source_doi': '10.1080/00207160.2014.946412',
        'arithmetic': 'fractions.Fraction; truncated power series modulo t**12',
        'path': '(t, 2*t)',
        'scope': 'One-step expansions and uniform local bounds, not the exact Q-order of every orbit.',
        'fields': {},
    }
    for label, case in CASES.items():
        calculation = analyze(case['coefficients'])
        if calculation['jacobian_determinant'][0] != 1 or calculation['sum_jacobians_determinant'][0] != 4:
            raise AssertionError('The required inverses are not justified at the origin.')
        values = {}
        for (name, vector), wanted in zip(calculation['outputs'].items(), case['expected']):
            actual = leading_term(vector)
            if actual != wanted:
                raise AssertionError(f'{label}/{name}: {actual} != {wanted}')
            values[name] = {
                'leading_degree': actual[0],
                'leading_coefficient': [str(actual[1]), str(actual[2])],
                'coefficients_through_degree_11': [[str(v) for v in component] for component in vector],
            }
            print(f'{label}/{name}: t^{actual[0]} * ({actual[1]}, {actual[2]})')
        report['fields'][label] = {'field': case['field'], 'methods': values,
                                 'jacobian_determinant_at_zero': '1',
                                 'sum_jacobians_determinant_at_zero': '4'}
        if label == 'coupled_quadratic':
            z = calculation['z']
            d = calculation['I_minus_K']
            if leading_term(z) != (3, Q(7, 2), Q(11, 2)):
                raise AssertionError('Incorrect leading coefficient of z.')
            expected_d = [[Q(-1), Q(1, 2)], [Q(1), Q(-1, 2)]]
            for i in range(2):
                for j in range(2):
                    if d[i][j][:3] != [Q(0), Q(0), expected_d[i][j]]:
                        raise AssertionError('Incorrect leading matrix of I-K.')
            report['proof_coefficients'] = {
                'z_degree_3': ['7/2', '11/2'],
                'I_minus_K_degree_2': [[str(v) for v in row] for row in expected_d],
            }
    report['status'] = 'All six leading coefficients and the analytic proof coefficients agree exactly.'
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'results' / 'lotfi_orders.json')
    args = parser.parse_args()
    report = verify()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(report['status'])


if __name__ == '__main__':
    main()

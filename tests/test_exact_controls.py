"""Regression tests for rational controls and exact truncated series."""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))
import gamma_controls as gamma
import verify_lotfi_orders as lotfi


class RationalControlTests(unittest.TestCase):
    def test_gamma_radius(self):
        self.assertTrue(gamma.in_rho(Q(3661, 12500)))
        self.assertFalse(gamma.in_rho(Q(3, 10)))
        self.assertFalse(gamma.in_rho(Q(-1, 10)))

    def test_upper_rounding(self):
        for x in (Q(0), Q(1, 3), Q(2, 7), Q(51, 100), Q(10, 3)):
            for digits in (0, 2, 24):
                value = gamma.up(x, digits)
                self.assertGreaterEqual(value, x)
                self.assertLess(value-x, Q(1, 10**digits))

    def test_zero_residual_is_fixed(self):
        state = gamma.cycle(Q(1, 10), Q(0))
        self.assertEqual(state['tau_next'], Q(1, 10))
        self.assertEqual(state['eps_next'], 0)
        for name in ('A', 'Y', 'B', 'W', 'C'):
            self.assertEqual(state[name], 0)

    def test_domain_rejection(self):
        with self.assertRaises(AssertionError):
            gamma.cycle(Q(3, 10), Q(1, 100))

    def test_first_cycle(self):
        state = gamma.initial_polynomials(Q(1547, 10000))
        self.assertEqual(state['V'], Q(292871921615, 10**12))
        self.assertLess(state['T1'], Q(1, 5))
        self.assertLess(state['P3'], Q(1, 75))

    def test_rational_transmitted_states(self):
        first = gamma.cycle(Q(1, 5), Q(1, 75))
        tau2, eps2 = gamma.up(first['tau_next']), gamma.up(first['eps_next'])
        self.assertEqual(tau2, Q(236509841095667281330063, 10**24))
        self.assertEqual(eps2, Q(426195191308547945528, 10**24))
        second = gamma.cycle(tau2, eps2)
        self.assertEqual(gamma.up(second['tau_next']), Q(238026135303284715478357, 10**24))
        self.assertEqual(gamma.up(second['eps_next']), Q(5818357024425655, 10**24))

    def test_tail(self):
        value = gamma.tail(Q(6,25), Q(238026136,10**9), Q(6,10**9), Q(1,2))
        self.assertEqual(value['m'], Q(97,361))
        self.assertEqual(value['M'], Q(31250,6859))
        self.assertLess(value['P_over_e'], Q(68,10**21))
        self.assertLess(value['tau_limit_bound'], Q(238026181,10**9))
        self.assertLess(value['stage_bound'], Q(238026159,10**9))


class SeriesTests(unittest.TestCase):
    def test_scalar_reciprocal(self):
        a = lotfi.constant(1)
        a[1], a[2] = Q(2), Q(-3)
        self.assertEqual(lotfi.multiply(a, lotfi.reciprocal(a)), lotfi.constant(1))

    def test_singular_reciprocal(self):
        with self.assertRaises(ZeroDivisionError):
            lotfi.reciprocal(lotfi.constant())

    def test_matrix_inverse(self):
        t = lotfi.constant(); t[1] = Q(1)
        matrix = [[lotfi.add(lotfi.constant(1), lotfi.scale(t,2)), t],
                  [lotfi.scale(t,4), lotfi.add(lotfi.constant(1), t)]]
        inverse = lotfi.matrix_inverse(matrix)
        self.assertEqual(lotfi.matrix_multiply(matrix, inverse), lotfi.identity())
        self.assertEqual(lotfi.matrix_multiply(inverse, matrix), lotfi.identity())

    def test_coupled_leading_terms(self):
        case = lotfi.CASES['coupled_quadratic']
        data = lotfi.analyze(case['coefficients'])
        self.assertEqual([lotfi.leading_term(v) for v in data['outputs'].values()], case['expected'])
        self.assertEqual(lotfi.leading_term(data['z']), (3,Q(7,2),Q(11,2)))

    def test_symmetric_leading_terms(self):
        case = lotfi.CASES['symmetric_quadratic']
        data = lotfi.analyze(case['coefficients'])
        self.assertEqual([lotfi.leading_term(v) for v in data['outputs'].values()], case['expected'])


if __name__ == '__main__':
    unittest.main()

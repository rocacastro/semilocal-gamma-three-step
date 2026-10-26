#!/usr/bin/env python3
"""Independent high-precision evaluation of the coupled Lotfi example.

This is a numerical consistency check, not the proof of the Taylor orders.
The exact proof coefficients are checked by verify_lotfi_orders.py.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import mpmath as mp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1]/'results'/'lotfi_numerical.json')
    args = parser.parse_args()
    rows = []
    with mp.workdps(200):
        def field(v):
            a, b = v
            return mp.matrix([a+a*b, b+a*a+a*b])
        def jacobian(v):
            a, b = v
            return mp.matrix([[1+b, a], [2*a+b, 1+a]])
        for exponent in (3, 5, 7):
            t = mp.mpf(10)**(-exponent)
            x = mp.matrix([t, 2*t]); jx = jacobian(x)
            y = x-mp.lu_solve(jx, field(x))
            jy = jacobian(y)
            z = x-2*mp.lu_solve(jx+jy, field(x))
            b = jx**-1*jy
            h = mp.mpf(7)/2*mp.eye(2)-4*b+mp.mpf(3)/2*b*b
            point = z
            quotients = {}
            for name, degree, coefficient in [('M6', 5, -mp.mpf(3)/4),
                                               ('M9', 7, mp.mpf(9)/8),
                                               ('M12', 9, -mp.mpf(27)/16)]:
                point = point-h*mp.lu_solve(jx, field(point))
                ratio = point/t**degree
                error = max(abs(ratio[0]-coefficient), abs(ratio[1]+coefficient))
                # Regression tolerance for direct multiprecision samples, not an interval bound.
                if error >= 1000*t:
                    raise ArithmeticError(f'Unexpected numerical discrepancy for {name}.')
                quotients[name] = {'scaled_step': [mp.nstr(v, 60) for v in ratio],
                                   'max_difference_from_leading_coefficient': mp.nstr(error, 30)}
            rows.append({'t': mp.nstr(t, 10), 'quotients': quotients})
    output = {'arithmetic': 'mpmath, 200 decimal digits, no intervals',
              'scope': 'Numerical consistency only; exact expansions are in lotfi_orders.json.',
              'samples': rows}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    print('All nine independent scaled-step checks passed (multiprecision, not intervals).')

if __name__ == '__main__':
    main()

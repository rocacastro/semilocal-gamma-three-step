#!/usr/bin/env python3
"""Additional checks for the TJM 2022 audit.

The bounds for Example 2.5 are verified with Fraction and integer powers.
The roots of the local functions are mpmath reproductions, not certified
intervals. The expansion of F(u)=u+u|u| is verified with SymPy.
Usage: python check_additional_2022.py --out additional_results.json
"""
from __future__ import annotations
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
import sympy as sp


def upper_power(x: Q, p: int, q: int, decimal_places: int = 20) -> Q:
    """Exact upper bound, with denominator 10**decimal_places, for x**(p/q)."""
    if x < 0 or p < 0 or q < 1:
        raise ValueError('Nonnegative x,p and positive q are required.')
    scale = 10**decimal_places
    target = x**p * scale**q
    lo, hi = 0, scale
    while Q(hi)**q < target:
        hi *= 2
    while hi - lo > 1:
        mid = (lo+hi)//2
        if Q(mid)**q < target:
            lo = mid
        else:
            hi = mid
    bound = Q(hi, scale)
    assert bound**q >= x**p
    return bound


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path(__file__).with_name('additional_2022.json'))
    args = parser.parse_args()
    data = {'exact': {}, 'mpmath_reproductions': {}}
    # Contraction margins certified by rational comparisons.
    for name, radius, displacement_limit, lip_lim in [
        ('small_ball',Q('0.238277'),Q('0.196'),Q('0.259')),
        ('large_ball',Q('2.54114'),None,Q('0.516'))]:
        u = 1+radius
        pow_35 = upper_power(u,3,5)
        lip = (Q(8,5)*pow_35+u/5)/8
        assert lip < lip_lim < 1
        info = {'radius':str(radius), 'rational_Lipschitz_bound':str(lip),
                'Lipschitz_less_than':str(lip_lim)}
        if displacement_limit is not None:
            displacement = (u*pow_35+u*u/10)/8
            assert displacement < displacement_limit < radius
            info.update(displacement_bound=str(displacement),displacement_less_than=str(displacement_limit))
        data['exact'][name]=info
    # The omega bound is refuted by constant functions of opposite signs.
    eps=Q(1,100)
    up_eps=upper_power(eps,3,5)
    lo_eps=up_eps-Q(1,10**20)
    up_dos=upper_power(2*eps,3,5)
    lhs_lo=Q(2,5)*lo_eps+eps/20
    rhs_up=(up_dos+2*eps)/5
    assert lhs_lo>rhs_up
    data['exact']['omega_difference']={'left_hand_side_lower_bound':str(lhs_lo),
                                        'right_hand_side_upper_bound':str(rhs_up)}
    # Branchwise Taylor expansion: y>0 and z<0 for small t>0.
    t=sp.symbols('t',positive=True)
    j=1+2*t
    y=sp.cancel(t-(t+t*t)/j)
    z=sp.cancel(y-5*(y+y*y)/j)
    assert sp.limit(z/t**2,t,0,dir='+')==-4
    step=sp.cancel(y-(9*(y+y*y)+(z-z*z))/(5*j))
    expansion=sp.series(step,t,0,5)
    assert sp.expand(expansion.removeO())==sp.Rational(32,5)*t**4
    data['exact']['C1_not_C2_operator']=str(expansion)
    # Numerical reproduction of radii, without claiming validated intervals.
    mp.mp.dps=70
    def omega(t): return (mp.mpf('1.5')*mp.sqrt(t)+t)/8
    def aux_integral(t): return 1+(mp.sqrt(t)+t/2)/8
    def g1(t): return (mp.sqrt(t)+t/2)/(8*(1-omega(t)))
    def g2(t):
        v=g1(t)
        return v+5*aux_integral(v*t)*v/(1-omega(t))
    def g3(t,coef):
        a,b=g1(t),g2(t)
        return b+(coef*aux_integral(a*t)*a+aux_integral(b*t)*b)/(5*(1-omega(t)))
    def bisect(f,lo,hi):
        lo,hi=mp.mpf(lo),mp.mpf(hi)
        assert f(lo)<0<f(hi)
        for _ in range(235):
            mid=(lo+hi)/2
            if f(mid)<0:lo=mid
            else:hi=mid
        return (lo+hi)/2
    vals={'r0':((mp.sqrt(137)-3)/4)**2,
          'r1':bisect(lambda t:g1(t)-1,'2','3'),
          'r2':bisect(lambda t:g2(t)-1,'.3','.6'),
          'r3_printed':bisect(lambda t:g3(t,1)-1,'.2','.4'),
          'r3_coefficient16':bisect(lambda t:g3(t,16)-1,'.15','.3')}
    for k,v in vals.items():data['mpmath_reproductions'][k]=mp.nstr(v,36)
    for radius in ['0.238277','2.54114']:
        u=1+mp.mpf(radius)
        data['mpmath_reproductions']['Lipschitz_'+radius]=mp.nstr(((mp.mpf(8)/5)*u**mp.mpf('.6')+u/5)/8,36)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Rational and symbolic checks: OK.')
    print('Numerical roots, not intervals:')
    for k,v in vals.items():print(k,mp.nstr(v,20))
    print('Output:',args.out)

if __name__=='__main__':main()

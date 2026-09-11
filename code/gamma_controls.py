#!/usr/bin/env python3
"""Exact controls for the gamma criterion; standard library only."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse, json, platform, sys


def up(x: F, digits: int = 24) -> F:
    scale = 10**digits
    return F(-((-x.numerator*scale)//x.denominator), scale)


def decimal(x: F, digits: int = 18, upper: bool = False) -> str:
    """Decimal representation using integers; upper=True rounds upwards."""
    assert x >= 0
    scale = 10**digits
    n = (-((-x.numerator*scale)//x.denominator)
         if upper else (x.numerator*scale)//x.denominator)
    return f'{n//scale}.{n%scale:0{digits}d}'


def in_rho(x: F, strict: bool = True) -> bool:
    """x < rho=1-1/sqrt(2) (or <=), decided exactly."""
    if not F(0) <= x < F(1):
        return False
    test = 2*(1-x)**2
    return test > 1 if strict else test >= 1


def phi(a: F) -> F:
    return a+5*a*a+5*a**3


def m(t: F) -> F:
    return 2-1/(1-t)**2


def M(t: F) -> F:
    return 2/(1-t)**3


def L(t: F, d: F) -> F:
    return d*(2*(1-t)-d)/((1-t)**2*(1-t-d)**2)


def E2(t: F, d: F) -> F:
    return d*d/((1-t)**2*(1-t-d))


def E3(t: F, d: F) -> F:
    return d**3/((1-t)**3*(1-t-d))


def initial_polynomials(a: F) -> dict[str, F]:
    """First-cycle bounds. These are polynomials in a with coefficients >=0.

    They can be evaluated algebraically outside the admissible interval. This
    does NOT assert that F is defined outside the gamma ball. The proof first
    establishes phi(alpha)<=rho for the actual parameter, then uses monotonicity
    of these polynomials to evaluate them at a rational upper bound for alpha.
    """
    Y=a*a+a**3
    B=5*Y
    ell=2*a+3*a*a
    e2=lambda t,d:(1+3*t)*d*d+d**3
    W=(ell*B+e2(a,B))/5
    C=Y+W
    P2=ell*W+e2(a,C)+e2(a,B)/5
    P3=ell*W+(2+6*a)*(W*W/2+W*B/5+F(2,25)*B*B)+C**3+B**3/5
    return dict(alpha=a,Y=Y,B=B,V=phi(a),ell=ell,W=W,C=C,
                T1=a+C,P2=P2,P3=P3)


def cycle(t: F, e: F) -> dict[str, F]:
    """Standard rational recurrence applied from n=1 after the polynomial initialization."""
    assert in_rho(t) and e>=0
    mn=m(t)
    assert mn>0
    A=e/mn; S=t+A
    assert in_rho(S)
    Y=E2(t,A); B=5*Y/mn; V=S+B
    assert in_rho(V)
    ell=L(t,A); W=(ell*B+E2(S,B))/(5*mn); C=B/5+W
    tn=S+C
    assert in_rho(tn)
    P2=ell*W+E2(S,C)+E2(S,B)/5
    P3=ell*W+M(S)*(W*W/2+W*B/5+F(2,25)*B*B)+E3(S,C)+E3(S,B)/5
    en=min(P2,P3)
    return dict(tau=t,eps=e,A=A,S=S,Y=Y,B=B,V=V,ell=ell,W=W,C=C,
                tau_next=tn,eps_next=en,P2=P2,P3=P3)


def tail(sigma: F, T: F, e: F, q: F) -> dict[str, F]:
    assert in_rho(sigma) and 0<q<1 and e>0
    mm,MM=m(sigma),M(sigma)
    a=e/mm; b=5*MM*a*a/(2*mm)
    w=(MM*a*b+MM*b*b/2)/(5*mm)
    c=b/5+w
    P=MM*a*w+MM*(c*c+b*b/5)/2
    bound=T+(a+c)/(1-q)
    stage=T+a+b
    assert P <= q*e
    assert max(bound,stage)<sigma
    return dict(sigma=sigma,T=T,e=e,q=q,m=mm,M=MM,a=a,b=b,w=w,c=c,
                P=P,P_over_e=P/e,tau_limit_bound=bound,stage_bound=stage)


def encode(d: dict[str,F]) -> dict:
    return {k: {'fraction':str(v), 'decimal_truncated':decimal(v),
                'decimal_upper':decimal(v,18,True)} for k,v in d.items()}


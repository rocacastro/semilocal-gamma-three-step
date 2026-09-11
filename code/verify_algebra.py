#!/usr/bin/env python3
"""Identities of the revision: SymPy, without numerical sampling tests."""
from pathlib import Path
import json
import sympy as s

def main():
    t,d,e,m,M,a=s.symbols('t d e m M a')
    checks={}
    k=lambda x:x*x/(1-x)
    kp=s.diff(k(t),t); kpp=s.diff(k(t),t,2)
    L=d*(2*(1-t)-d)/((1-t)**2*(1-t-d)**2)
    E2=d*d/((1-t)**2*(1-t-d))
    E3=d**3/((1-t)**3*(1-t-d))
    checks['L_integrated']=s.cancel(L-(kp.subs(t,t+d)-kp))==0
    checks['E2_integrated']=s.cancel(E2-(k(t+d)-k(t)-kp*d))==0
    checks['E3_integrated']=s.cancel(E3-(E2-kpp*d*d/2))==0
    ca=e/m; cb=5*M*ca*ca/(2*m)
    cw=(M*ca*cb+M*cb*cb/2)/(5*m);cc=cb/5+cw
    P=M*ca*cw+M*(cc*cc+cb*cb/5)/2
    expected=(5*M**3*e**4/(4*m**6)+7*M**4*e**5/(8*m**8)
        +7*M**5*e**6/(16*m**10)+5*M**6*e**7/(16*m**12)
        +25*M**7*e**8/(128*m**14))
    checks['tail_polynomial_P']=s.expand(P-expected)==0
    g=lambda x:a-x+x*x/(1-x)
    v=a+5*g(a)
    t1=s.factor(a+(9*g(a)+g(v))/5)
    checks['first_majorant_iterate']=s.cancel(t1-a*(1-2*a)*(1-2*a*a)/((1-a)*(1-2*a-4*a*a)))==0
    exact=2*a**5*(5-4*a-8*a*a+8*a**3)/((1-a)*(1-2*a-4*a*a)*(1-4*a+6*a**3-4*a**4))
    checks['rational_majorant_residual']=s.cancel(g(t1)-exact)==0
    checks['residual_coefficient_10']=s.limit(exact/a**5,a,0)==10
    h=s.symbols('h');x=s.Matrix([h,2*h])
    F=lambda x:s.Matrix([x[0]+x[1]**2,x[1]+x[0]**2])
    J=s.Matrix([[1,4*h],[2*h,1]])
    inv=J.inv();y=x-inv*F(x);z=y-5*inv*F(y)
    xn=y-inv*(9*F(y)+F(z))/5
    for i,coeff in enumerate((28,-56)):
        checks[f'fourth_degree_coefficient_component_{i+1}']=s.limit(xn[i]/h**4,h,0)==coeff
    checks['all']=all(checks.values())
    if not checks['all']:raise ArithmeticError(checks)
    out=Path(__file__).resolve().parent.parent/'results'
    out.mkdir(exist_ok=True)
    (out/'algebra_verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    for name,value in checks.items():print(name,value)

if __name__=='__main__':main()

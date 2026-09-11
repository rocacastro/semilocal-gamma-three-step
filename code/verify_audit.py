#!/usr/bin/env python3
"""Reproducible audit of Cordero et al. (JCAM 2015) and Singh et al. (AMC 2016).

Sources: DOI 10.1016/j.cam.2014.06.008 and 10.1016/j.amc.2015.11.062.
Dependencies: sympy and mpmath. Python >= 3.10.
Usage: python verify_audit.py --dps 2000 --out results

Logical distinction:
- The Fraction/SymPy checks are exact.
- The mpmath tables are numerical reproductions, NOT interval certificates.
- A false inequality does not by itself refute convergence of the method.
- The 2016 local counterexample specifically refutes uniqueness in the printed
  radius. The corrected radius includes the power 1/q.
- The vector expansion refutes a UNIFORM local bound O(||e||^5); it does not
  assert that all trajectories have identical Q-order.
"""
from __future__ import annotations
import argparse
import csv
import json
from fractions import Fraction as Q
from pathlib import Path
try:
    import sympy as sp
    import mpmath as mp
except ImportError as exc:
    raise SystemExit('Install: python -m pip install sympy mpmath') from exc


def exact_checks() -> dict:
    t = sp.symbols('t', real=True)
    def trunc(v):
        return v.applyfunc(lambda c: sp.series(c, t, 0, 6).removeO().expand())
    def field(v):
        u, w = v
        return sp.Matrix([u+w*w, w+u*u])
    x = sp.Matrix([t, 2*t])
    inv = sp.Matrix([[1,4*t],[2*t,1]]).inv()
    y = trunc(x-inv*field(x))
    z = trunc(y-5*inv*trunc(field(y)))
    xp = trunc(z-inv*(-16*trunc(field(y))+trunc(field(z)))/5)
    xe = trunc(y-inv*(9*trunc(field(y))+trunc(field(z)))/5)
    expected = sp.Matrix([28*t**4+224*t**5, -56*t**4+70*t**5])
    assert xp == xe == expected

    # The 2015 h retains exactly the factor 5 from the second stage.
    a = sp.symbols('a', positive=True)
    h = a/2+a*a/2+5*a**3/8
    assert sp.expand(a*(4+(1+5*a/2)**2)/10-h) == 0
    F = lambda u:u-1-a*u*u/2
    ys = sp.Integer(1); zs = sp.expand(ys-5*F(ys))
    xs = sp.expand(ys-(9*F(ys)+F(zs))/5)
    assert xs == 1+h  # in the concave example, bound (5) from 2015 is exact

    # False INTERMEDIATE inequality in (6), 2015.
    # F(u,v)=(u,100v), x0=(0,1), ||Gamma||=beta=eta=1, K=1/10.
    aa = Q(1,10)
    hh = aa/2+aa*aa/2+Q(5,8)*aa**3
    assert Q(100) > 1+hh

    # 2016: C-infinity example refuting (1+t(c0))^q d0 in (2.17).
    # F(u)=u-1-u^2/2000, D=(-1.01,1.01), beta=eta=1,
    # q=1/2, K=c0=1/700. The semilocal criterion and the domain ARE satisfied.
    c=Q(1,700); w=Q(1,1050)
    # t(c)=w+w^2+w(1+5w)/sqrt(5250); sqrt(5250)>72.
    assert Q(72)**2 < 5250
    Tupper=w+w*w+w*(1+5*w)/72
    T=Q(967,10**6)
    assert Tupper<T
    x1=Q(1)+Q(1,2000)+Q(1,2*10**6)+Q(5,8*10**9)
    assert x1 == Q('1.000500500625')
    assert x1 > 1+T/2  # sqrt(1+t) <= 1+t/2 < 1+T/2
    assert Q(1,1000)**2*Q(202,100) < c*c  # Holder on D
    assert T < Q(32,1000)**2
    supper=w+(1+c)*T+w*T*Q(32,1000)
    rupper=1/(1-c*(1+T/2))
    assert rupper**3*supper < 1  # r^(1+q)s^q<1, q=1/2
    Rupper=Q(5)*c/Q(3,2)+(1+T)/(1-rupper*supper)
    assert Rupper < Q(101,100)  # semilocal ball contained in D

    # LOCAL Theorem 3.1 from 2016: false uniqueness radius due to omission of 1/q.
    # G(u)=u-4u^2, x*=0, D0=(-.3,.3), q=.5, K0=9/2, K1=25/4, R=7/25.
    K0=Q(9,2);K1=Q(25,4);D=Q(3,10);R=Q(7,25)
    assert Q(8)**2*D < K0*K0
    assert Q(8)**2*(2*D) < K1*K1
    assert Q(1,4) < R < Q(3,2)/K0
    assert R<D
    assert Q(1,4)-4*Q(1,4)**2 == 0
    # r1=( (1+q)/(K1+(1+q)K0) )^(1/q), and rho<r1.
    r1=(Q(3,2)/(K1+Q(3,2)*K0))**2
    assert r1<R
    corrected=(Q(3,2)/K0)**2
    assert corrected == Q(1,9) and corrected < Q(1,4)

    # Local Example 3.2: 2(sqrt(3)-1)>sqrt(2); equivalent to 196>192.
    assert 196>192

    return {
        'vector_expansion_through_degree_5':str(expected.T),
        '2015_h_retains_factor_5':True,
        '2015_equation6_false':{'lhs':100,'rhs':str(1+hh)},
        '2016_equation217_false':{'exact_x1':str(x1),'printed_bound_less_than':str(1+T/2),
                               'certified_R_less_than':'1.01'},
        '2016_local_uniqueness_refuted':{'roots':['0','1/4'],'R':'7/25','K0':'9/2','K1':'25/4',
                                      'q':'1/2','correct_limit':'1/9'}}



def functions(c,q):
    t=c/(1+q)+c*c/(1+q)**2+mp.power(5,q-1)*c**(1+q)/(1+q)**(1+q)+mp.power(5,q)*c**(2+q)/(1+q)**(2+q)
    r=1/(1-c*(1+t)**q)
    s=c/(1+q)+(c+1)*t+c*t**(1+q)/(1+q)
    return r,s,t


def monotone_root(fun,lo,hi,steps=240):
    for _ in range(steps):
        mid=(lo+hi)/2
        if fun(mid)>0:hi=mid
        else:lo=mid
    return (lo+hi)/2


def threshold(q):
    def fun(c):
        r,s,t=functions(c,q)
        den=1-c*(1+t)**q
        return mp.mpf(1) if den<=0 else q*mp.log(s)-(1+q)*mp.log(den)
    return monotone_root(fun,mp.mpf('1e-80'),mp.mpf('.7'))


def local_radii(k0,k1,q):
    def gam(t):
        D=1-k0*t**q
        g1=k1*t**q/((1+q)*D)
        g2=(1+5*(1+k0*g1**q*t**q)/D)*g1
        g3=g2+(16*(1+k0*g1**q*t**q)*g1+(1+k0*g2**q*t**q)*g2)/(5*D)
        return g1,g2,g3
    r1=((1+q)/(k1+(1+q)*k0))**(1/q)
    r2=monotone_root(lambda t:gam(t)[1]-1,mp.mpf(0),r1)
    rho=monotone_root(lambda t:gam(t)[2]-1,mp.mpf(0),r2)
    return r1,r2,rho


def singh_numerics() -> dict:
    with mp.workdps(80):
        fmt=lambda x:mp.nstr(x,30)
        qs=[mp.mpf(v) for v in ['.85','1','.7','.5']]
        result={'thresholds':{str(q):fmt(threshold(q)) for q in qs},'semilocal':[],'local':[]}
        # Corrected constants following from the operators in the examples;
        # the printed formulas for eta contain misprints identified in the report.
        for name,q,K in [('2.1, q=.85',mp.mpf('.85'),mp.mpf('1.85')/24),
                          ('2.1, q=1',mp.mpf(1),mp.mpf(1)/12),
                          ('2.2, q=.7',mp.mpf('.7'),mp.mpf('1.7')*mp.log(2)/4),
                          ('2.2, q=1',mp.mpf(1),mp.log(2)/2)]:
            beta=1/(1-K)
            eta=beta/24 if name.startswith('2.1') else beta*mp.log(2)/4
            c=K*beta*eta**q;r,s,t=functions(c,q)
            R=5*c/(1+q)+(1+t)/(1-r*s)
            u=eta*((1+q)/c-R**q)**(1/q)
            Rwrong=5*c/(1+q)+(1+t)**q/(1-r*s)
            uwrong=eta*((1+q)/c-(Rwrong*eta)**q)**(1/q)
            result['semilocal'].append({'example':name,'beta':fmt(beta),'eta':fmt(eta),'c':fmt(c),
                                         'existence_radius_theorem_formula':fmt(R*eta),
                                         'uniqueness_radius_theorem_formula':fmt(u),
                                         'existence_radius_with_incorrect_power':fmt(Rwrong*eta),
                                         'uniqueness_radius_with_double_scaling':fmt(uwrong)})
        for name,k0,k1,q in [('3.1','96.6628','96.6628','.4'),('3.2','1','1','.5'),
                            ('3.3',mp.e-1,mp.e,'.7'),('3.4','7.5','15','.5')]:
            k0,k1,q=map(mp.mpf,[k0,k1,q])
            result['local'].append({'example':name,'q':str(q),
                'r1_r2_rho':[fmt(v) for v in local_radii(k0,k1,q)],
                'rho_with_q_05':fmt(local_radii(k0,k1,mp.mpf('.5'))[2])})
        # Detailed value of counterexample (2.17), supported above by exact bounds.
        r,s,t=functions(mp.mpf(1)/700,mp.mpf('.5'))
        result['counterexample_217']={'sqrt_1_plus_t':fmt(mp.sqrt(1+t)),
                                   'radius_R':fmt(mp.mpf(1)/210+(1+t)/(1-r*s))}
        return result


def reproduce_2015(dps:int,out:Path) -> dict:
    with mp.workdps(dps):
        n=8
        nodes,weights=mp.gauss_quadrature(n,'legendre')
        nodes=[(nodes[j]+1)/2 for j in range(n)]
        weights=[weights[j]/2 for j in range(n)]
        A=mp.matrix(n)
        for i in range(n):
            for j in range(n):
                A[i,j]=weights[j]*min(nodes[i],nodes[j])*(1-max(nodes[i],nodes[j]))
        I=mp.eye(n);one=mp.matrix([1]*n)
        def vn(v):return max(abs(t) for t in v)
        def mn(B):return max(sum(abs(B[i,j]) for j in range(B.cols)) for i in range(B.rows))
        def F(x):return x-one-A*mp.matrix([v*v for v in x])
        def J(x):return I-2*A*mp.diag(list(x))
        fmt=lambda x:mp.nstr(x,35)
        x0=mp.matrix([mp.mpf('1.6')]*n)
        K=2*mn(A);beta=mn(J(x0)**-1);eta=vn(mp.lu_solve(J(x0),F(x0)));a=K*beta*eta
        f,g,h=functions(a,mp.mpf(1))
        R=mp.mpf('2.5')*a+(1+h)/(1-f*g)
        x=x0.copy();hist=[];ds=[]
        for k in range(1,7):
            j=J(x)
            y=x-mp.lu_solve(j,F(x));z=y-5*mp.lu_solve(j,F(y))
            xx=y-mp.lu_solve(j,9*F(y)+F(z))/5
            d=vn(xx-x);res=vn(F(xx));ds.append(d)
            p=mp.log(ds[-1]/ds[-2])/mp.log(ds[-2]/ds[-3]) if k>=3 else None
            hist.append({'n':k,'residual':fmt(res),'increment':fmt(d),'ACOC':fmt(p) if p is not None else '',
                         'residual_at_working_precision_floor':res < mp.power(10, -dps+20)})
            x=xx
        with (out/'cordero2015_iterations.csv').open('w',newline='',encoding='utf-8') as fp:
            wr=csv.DictWriter(fp,fieldnames=hist[0].keys());wr.writeheader();wr.writerows(hist)
        return {'dps':dps,'N':8,'x0':'1.6 * (1,...,1)','K':fmt(K),'beta':fmt(beta),'eta':fmt(eta),'a0':fmt(a),
                'existence_radius':fmt(R*eta),'uniqueness_radius':fmt(2/(K*beta)-R*eta),
                'history':hist,'note':'Residuals close to 10^(-dps) are NOT certified error bounds.'}


def indices() -> dict:
    with mp.workdps(60):
        entries=[('Newton',2,8,368),('Jarratt',4,4,728),('WKG',6,3,1096),('M5_published_label',5,4,640),('M5_generic_order',4,4,640)]
        return {name:{'cost_per_step':C,'total_cost':k*C,'p_per_step_cost':mp.nstr(mp.mpf(p)**(1/mp.mpf(C)),18),
                      'p_per_total_cost':mp.nstr(mp.mpf(p)**(1/mp.mpf(k*C)),18)} for name,p,k,C in entries}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps',type=int,default=2000)
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results')
    args=parser.parse_args()
    if args.dps<1600:parser.error('Use at least 1600 digits for the six updates.')
    args.out.mkdir(parents=True,exist_ok=True)
    result={'scope':'Exact checks and numerical reproductions, distinguished. Does not establish priority.',
            'exact':exact_checks(),'singh2016_numerical':singh_numerics(),
            'cordero2015_numerical':reproduce_2015(args.dps,args.out),'indices_2015':indices()}
    path=args.out/'audit_results.json'
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print('All EXACT assertions were verified.')
    print('Numerical results (not intervals):',path)
    print('Last ACOC values of the 2015 reproduction:')
    for row in result['cordero2015_numerical']['history']:
        if row['ACOC']:print(row['n'],row['ACOC'])
    print('The current radial comparison is run in verify_gamma_criterion.py.')

if __name__=='__main__':main()

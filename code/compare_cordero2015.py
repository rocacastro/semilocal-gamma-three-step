#!/usr/bin/env python3
"""S3: rational comparison with Cordero et al. (2015).
The bounds in this table use the standard rational controls tau_n,
started at (0,alpha). They are not confused with the hatted sequences.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json
from gamma_controls import cycle,up,in_rho,tail as full_tail,decimal

def rat(s):return Q(str(s))
def dec(x,places=26):return decimal(x,places,True)
def ceil_grid(x,places=24):return up(x,places)
def below_gamma_radius(t):return in_rho(t)
def step(t,e):
    c=cycle(t,e)
    return {**c,'t_next':c['tau_next'],'e_next':c['eps_next']}
def tail(S,e):
    from gamma_controls import m,M
    mm,MM=m(S),M(S)
    a=e/mm;b=5*MM*a*a/(2*mm);w=(MM*a*b+MM*b*b/2)/(5*mm);c=b/5+w
    p=MM*a*w+MM*(c*c+b*b/5)/2
    return dict(a=a,b=b,c=c,p=p)
def h15(a): return a/2+a*a/2+5*a**3/8

def f15(a): return 1/(1-a*(1+h15(a)))

def g15(a): return a/2+(1+a)*h15(a)+a*h15(a)**2/2

def R15(a): return Q(5,2)*a+(1+h15(a))/(1-f15(a)*g15(a))

def comparison_certificate():
    lo=rat('0.29313992158369'); hi=rat('0.29313992158370')
    assert f15(lo)**2*g15(lo)<1<f15(hi)**2*g15(hi)
    # Polynomial witness F(u,v)=(u+v^2-3/20, v+u^2-3/40).
    a=Q(3,10)
    assert f15(a)**2*g15(a)>1
    return {'root_2015_interval':[str(lo),str(hi)],
            'witness_alpha':'3/20','witness_a_2015':str(a),
            'witness_f2g':dec(f15(a)**2*g15(a))}

def certified_localization(alpha, target=Q(1,10**24)):
    """Finite radial evolution + rigorous geometric tail for table values."""
    t=Q(0);e=alpha;states=[]
    for n in range(10):
        row=step(t,e);states.append(row)
        t=ceil_grid(row['t_next'],40); e=ceil_grid(row['e_next'],40)
        if e<target:break
    # A ball beyond t; tail test is exact. No assertion that this radius is optimal.
    S=ceil_grid(t+Q(1,10**7),8)
    assert below_gamma_radius(S)
    tt=tail(S,e);q=Q(1,2)
    assert tt['p']<q*e
    bound=t+(tt['a']+tt['c'])/(1-q)
    assert t+max((tt['a']+tt['c'])/(1-q),tt['a']+tt['b'])<S
    return bound, states

def table_certificate():
    # For the polynomial class ||F''||=2, all constants of the 2015 theorem
    # are exact, unlike a generic bound extracted from gamma.
    rows=[]
    for text in ['0.01','0.02','0.05','0.1','0.14','0.15']:
        alpha=rat(text)
        bnd,_=certified_localization(alpha)
        a=2*alpha; applicable=f15(a)**2*g15(a)<1
        row={'alpha':text,'radial_location_upper':dec(bnd),'2015_applies':applicable}
        if applicable:
            loc15=(1+h15(a))*alpha/(1-f15(a)*g15(a))
            ball15=R15(a)*alpha
            a1=a*f15(a)**2*g15(a)
            d1=f15(a)*g15(a)*alpha
            prefix=alpha*(1+h15(a))+d1*(1+h15(a1))
            assert bnd<prefix<loc15<ball15
            row.update({'2015_solution_location':dec(loc15),'2015_stage_radius':dec(ball15),
                        '2015_two_term_prefix':dec(prefix)})
        rows.append(row)
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent.parent/'results')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    report={'arithmetic':'fractions.Fraction, all comparisons exact','comparison':comparison_certificate(),'localization_table':table_certificate()}
    for row,printed in zip(report['localization_table'][2:],['0.053036119','0.116956871','0.195079730','0.230706676']):
        bound,_=certified_localization(rat(row['alpha']));assert bound<rat(printed)
        row['table_upper_bound']=printed
    (a.out/'comparison_2015.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print('All rational comparisons with 2015: OK.')
    for row in report['localization_table']:print(row)
if __name__=='__main__':main()

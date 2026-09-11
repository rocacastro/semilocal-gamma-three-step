#!/usr/bin/env python3
"""S3: exact rational verification of alpha <= 1547/10000.
Run: python verify_gamma_criterion.py --out results
No acceptance comparison uses floating point.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse,json,platform,sys
from gamma_controls import initial_polynomials,cycle,tail,up,decimal,in_rho,encode

def main():
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parent.parent/'results')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    alpha=F(1547,10000);sa=F(3661,12500);st=F(6,25)
    assert (3-alpha)**2-8==F(9573209,100000000)>0
    assert 2*(1-sa)**2-1==F(2921,78125000)>0
    c0=initial_polynomials(alpha)
    assert c0['V']==F(292871921615,10**12)<sa
    assert c0['T1']<F(198986457,10**9)<F(1,5)
    assert c0['P3']<F(13187894,10**9)<F(1,75)
    assert in_rho(sa) and in_rho(st)
    t,e=F(1,5),F(1,75);states=[]
    for n in (1,2):
        c=cycle(t,e)
        assert max(c['V'],c['tau_next'])<F(253,1000)<sa
        tn,en=up(c['tau_next']),up(c['eps_next'])
        states.append({'n':n,'values':encode(c),'upper_next_tau':str(tn),'upper_next_epsilon':str(en)})
        t,e=tn,en
    T=F(238026136,10**9);Et=F(6,10**9)
    assert t<T and e<Et
    tt=tail(st,T,Et,F(1,2))
    assert tt['P_over_e']<F(68,10**21)<F(1,2)
    assert tt['tau_limit_bound']<F(238026181,10**9)<st
    assert tt['stage_bound']<F(238026159,10**9)<st
    # Also verifies all presentation roundings in the table.
    bounds=[['0.252029367','0.236509842'],['0.238096744','0.238026136']]
    for row,b in zip(states,bounds):
        assert F(row['values']['V']['fraction'])<F(b[0])
        assert F(row['values']['tau_next']['fraction'])<F(b[1])
    result={'alpha_c':str(alpha),'sigma_all':str(sa),'sigma_tail':str(st),
      'initial_polynomials':encode(c0),'upper_initialization':{'tau1':'1/5','epsilon1':'1/75'},
      'cycles':states,'tail':encode(tt),'arithmetic':'fractions.Fraction; all inequalities exact',
      'status':'ALL EXACT COMPARISONS ARE TRUE','python':sys.version,'platform':platform.platform()}
    (args.out/'criterion_alpha1547.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    lines=[result['status'],'alpha_c = 1547/10000 = 0.1547',
       'Vhat_0 = '+decimal(c0['V'],24), 'tauhat_1 < '+decimal(c0['T1'],24,True),
       'epsilonhat_1 < '+decimal(c0['P3'],24,True)]
    for row in states:
        lines.append('n='+str(row['n'])+': '+ '; '.join(k+' <= '+row['values'][k]['decimal_upper'] for k in ['tau','eps','V','tau_next','eps_next']))
        lines.append('next data tau='+row['upper_next_tau']+'; epsilon='+row['upper_next_epsilon'])
    lines.extend(k+' <= '+decimal(tt[k],30,True) for k in ['P_over_e','tau_limit_bound','stage_bound'])
    text='\n'.join(lines)+'\n';(args.out/'criterion_alpha1547.txt').write_text(text);print(text)
if __name__=='__main__':main()

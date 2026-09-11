#!/usr/bin/env python3
"""S5: five cycles with one LU per cycle, Nystrom/Gauss--Legendre.

Python >= 3.10; mpmath 1.3.0. All parameters listed in Supplementary Section S5 are exact.
python reproduce_experiments.py --problem all --N 30 --out results/experiments
python reproduce_experiments.py --problem chandrasekhar --N 20 40 60 --out results/experiments
python reproduce_experiments.py --problem bratu --N 20 40 60 --out results/experiments
python reproduce_experiments.py --summarize --out results/experiments

Stopping uses a fixed number of cycles. Iterations use 1000 digits; profiles 100.
Residuals and grid maxima are NOT quadrature error bounds or results of
interval arithmetic. ACOC values near saturation are omitted.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,platform,sys,time
from pathlib import Path
try:
    import mpmath as mp
except ImportError as exc:
    raise SystemExit('Install: python -m pip install mpmath==1.3.0') from exc
PARAMS={'chandrasekhar':{'omega0':'1/2','x0':'1'},'polynomial':{'epsilon':'1','x0':'0'},
 'bratu':{'lambda':'2','x0':'0'},'nonseparable':{'rho_H':'1','lambda':'1/100','x0':'0'},
 'reciprocal':{'mu':'1/10','x0':'1'}}

def text(x,d=25):
    return '' if x is None else mp.nstr(x,d,min_fixed=0,max_fixed=0)

def csv_out(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def norm(v):return max((abs(z) for z in v),default=mp.mpf(0))

def problem(name,N):
    t,w=mp.gauss_quadrature(N,'legendre');t=[(z+1)/2 for z in t];w=[z/2 for z in w]
    if name in ('chandrasekhar','reciprocal'):
        K=mp.matrix([[w[j]*t[i]/(t[i]+t[j]) for j in range(N)] for i in range(N)])
    elif name=='bratu':
        K=mp.matrix([[w[j]*(t[j]*(1-t[i]) if t[j]<=t[i] else t[i]*(1-t[j])) for j in range(N)] for i in range(N)])
    elif name=='nonseparable':
        K=mp.matrix([[w[j]*mp.exp(t[i]*t[j]) for j in range(N)] for i in range(N)])
    else: coeff=[w[j]*t[j]**4 for j in range(N)]
    f0=mp.matrix([mp.sin(s/3) for s in t])
    def F(x):
        if name=='chandrasekhar':
            kx=K*x;return mp.matrix([x[i]-1-x[i]*kx[i]/4 for i in range(N)])
        if name=='reciprocal':
            if any(z<=0 for z in x):raise ArithmeticError('Outside the positive domain')
            v=K*mp.matrix([1/z for z in x]);return mp.matrix([x[i]-1+v[i]/10 for i in range(N)])
        if name=='bratu':return x-2*(K*mp.matrix([mp.exp(z) for z in x]))
        if name=='nonseparable':return x-f0-(K*mp.matrix([mp.expm1(z)-z for z in x]))/100
        integral=mp.fsum(coeff[j]*(x[j]**3+x[j]**2/2) for j in range(N))
        return mp.matrix([x[i]-f0[i]-t[i]*integral for i in range(N)])
    def J(x):
        A=mp.eye(N)
        if name=='chandrasekhar':
            kx=K*x
            for i in range(N):
                for j in range(N):A[i,j]-=x[i]*K[i,j]/4
                A[i,i]-=kx[i]/4
        elif name in ('reciprocal','bratu','nonseparable'):
            if name=='reciprocal':
                if any(z<=0 for z in x):raise ArithmeticError('Outside the positive domain')
                r=[1/(10*z*z) for z in x]
            elif name=='bratu':r=[2*mp.exp(z) for z in x]
            else:r=[mp.expm1(z)/100 for z in x]
            for i in range(N):
                for j in range(N):A[i,j]-=K[i,j]*r[j]
        else:
            r=[coeff[j]*(3*x[j]**2+x[j]) for j in range(N)]
            for i in range(N):
                for j in range(N):A[i,j]-=t[i]*r[j]
        return A
    x0=mp.matrix([1 if name in ('chandrasekhar','reciprocal') else 0]*N)
    return t,w,F,J,x0

def solve(fac,b):
    """P A=L U; substitutions only, without refactorization."""
    P,L,U=fac;pb=P*b;n=b.rows;y=mp.matrix(n,1);x=mp.matrix(n,1)
    for i in range(n):y[i]=(pb[i]-mp.fsum(L[i,j]*y[j] for j in range(i)))/L[i,i]
    for i in range(n-1,-1,-1):
        if U[i,i]==0:raise ArithmeticError('Zero pivot')
        x[i]=(y[i]-mp.fsum(U[i,j]*x[j] for j in range(i+1,n)))/U[i,i]
    return x

def run(name,N,dps,cycles,out):
    if N<2 or dps<80 or cycles<3:raise ValueError('N>=2, dps>=80 and cycles>=3 are required')
    dest=out/f'{name}_N{N}_dps{dps}';dest.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    with mp.workdps(dps):
        t,w,F,J,x=problem(name,N)
        moments=[abs(mp.fsum(w[j]*t[j]**k for j in range(N))-mp.mpf(1)/(k+1)) for k in (0,1,2,2*N-1)]
        if max(moments)>mp.power(10,-dps+20):raise ArithmeticError('Quadrature moment error')
        steps=[];history=[];vectors=[];linear_res=[]
        def save(n,xx,step=None,stages=None):
            r=norm(F(xx));floor=max(mp.mpf(1),norm(xx))*mp.power(10,-dps+30)
            p=None;status='fewer_than_three_increments'
            if n>=3:
                ds=steps[-3:]
                if min(ds)<=floor:status='zero_increment_or_near_precision_limit'
                else:
                    den=mp.log(ds[1]/ds[0])
                    if den==0:status='zero_denominator'
                    else:p=mp.log(ds[2]/ds[1])/den;status='usable'
            history.append(dict(n=n,residual_inf=text(r),increment_inf=text(step),acoc=text(p,16),
                                acoc_status=status,residual_near_precision=bool(r<=floor)))
            item={'n':n,'x':[text(z,dps+8) for z in xx]}
            if stages:item.update({k:[text(z,dps+8) for z in v] for k,v in stages.items()})
            vectors.append(item)
        save(0,x)
        for n in range(cycles):
            Fx=F(x);A=J(x);fac=mp.lu(A)
            p=solve(fac,Fx);y=x-p;Fy=F(y);q=solve(fac,Fy);z=y-5*q
            Fz=F(z);r=solve(fac,9*Fy+Fz);xn=y-r/5
            linear_res.append(text(max(norm(A*p-Fx),norm(A*q-Fy),norm(A*r-(9*Fy+Fz)))))
            step=norm(xn-x);steps.append(step);x=xn;save(n+1,x,step,{'previous_y':y,'previous_z':z})
        record={'problem':name,'N':N,'precision_digits':dps,'cycles':cycles,'parameters':PARAMS[name],
                'quadrature':'Gauss--Legendre [0,1], constructed at the working precision','norm':'maximum',
                'stopping_rule':'fixed number of cycles','linear_solver':'mpmath.lu, pivoting, one factorization per cycle',
                'rounding_model':'multiprecision floating point, without intervals','history':history,
                'quadrature_moment_errors':[text(e) for e in moments],'linear_equation_residuals':linear_res,
                'nodes':[text(z,dps+8) for z in t],'weights':[text(z,dps+8) for z in w],
                'solution':[text(z,dps+8) for z in x],'python':sys.version,'mpmath':mp.__version__,
                'backend':mp.libmp.BACKEND,'platform':platform.platform(),'machine':platform.machine(),
                'elapsed_seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (dest/'result.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        (dest/'complete_iterates.json').write_text(json.dumps(vectors,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        csv_out(dest/'history.csv',history)
        csv_out(dest/'nodes_weights_solution.csv',[dict(j=j+1,t=text(t[j],dps+8),w=text(w[j],dps+8),x=text(x[j],dps+8)) for j in range(N)])
    print(name,N,'r4=',history[4]['residual_inf'] if cycles>=4 else '-',
          'p5=',history[5]['acoc'] if cycles>=5 else '-',flush=True)

def profile(data,grid):
    t=list(map(mp.mpf,data['nodes']));w=list(map(mp.mpf,data['weights']));x=list(map(mp.mpf,data['solution']))
    if data['problem']=='chandrasekhar':
        return [mp.mpf(1) if s==0 else 1/(1-mp.fsum(wj*s*xj/(s+tj) for tj,wj,xj in zip(t,w,x))/4) for s in grid]
    if data['problem']=='bratu':
        ex=list(map(mp.exp,x))
        return [2*mp.fsum(wj*(tj*(1-s) if tj<=s else s*(1-tj))*ej for tj,wj,ej in zip(t,w,ex)) for s in grid]
    raise ValueError('Only Chandrasekhar and Bratu have a reconstruction here')

def summarize(out):
    data={}
    for p in sorted(out.glob('*/result.json')):
        d=json.loads(p.read_text());data[d['problem'],d['N']]=d
    baseline=[]
    for name in PARAMS:
        d=data[name,30]
        for h in d['history']:baseline.append({'problem':name,'N':30,'dps':d['precision_digits'],**h})
    csv_out(out/'five_problems_N30.csv',baseline)
    rows=[];profiles=[]
    with mp.workdps(100):
        grid=[mp.mpf(i)/400 for i in range(401)]
        a=mp.findroot(lambda a:2*a*a/mp.cosh(a/2)**2-2,(mp.mpf(1),mp.mpf('1.5')))
        assert 0<a<2 and abs(2*a*a/mp.cosh(a/2)**2-2)<mp.mpf('1e-95')
        exact=[2*mp.log(mp.cosh(a/2)/mp.cosh(a*(s-mp.mpf(1)/2))) for s in grid]
        for name in ('chandrasekhar','bratu'):
            ref=profile(data[name,60],grid)
            for N in (20,30,40,60):
                d=data[name,N];pr=profile(d,grid)
                e60=max(abs(u-v) for u,v in zip(pr,ref))
                ex=max(abs(u-v) for u,v in zip(pr,exact)) if name=='bratu' else None
                rows.append({'problem':name,'N':N,'iteration_dps':d['precision_digits'],'profile_dps':100,
                  'acoc_5':d['history'][5]['acoc'],'residual_4':d['history'][4]['residual_inf'],
                  'difference_vs_N60':text(e60,22),'difference_vs_exact_Bratu':text(ex,22)})
                profiles += [{'problem':name,'N':N,'i':i,'s':text(s,10),'u_N':text(u,100),
                              'u_exact_Bratu':text(exact[i],100) if name=='bratu' else ''}
                             for i,(s,u) in enumerate(zip(grid,pr))]
        csv_out(out/'refinement.csv',rows);csv_out(out/'profiles_grid401.csv',profiles)
        (out/'Bratu_reference.json').write_text(json.dumps({'lambda':'2','a':text(a,100),
          'formula':'u(s)=2 log(cosh(a/2)/cosh(a(s-1/2)))','equation':'2a^2/cosh(a/2)^2=2',
          'branch':'smallest positive root, 0<a<2','precision':100},indent=2)+'\n')
    print('Summary and profiles generated in',out,flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--problem',choices=['all',*PARAMS],default='all');p.add_argument('--N',type=int,nargs='+',default=[30])
    p.add_argument('--dps',type=int,default=1000);p.add_argument('--cycles',type=int,default=5)
    p.add_argument('--out',type=Path,default=Path('results/experiments'));p.add_argument('--summarize',action='store_true')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if a.summarize:summarize(a.out);return
    for name in PARAMS if a.problem=='all' else [a.problem]:
        for N in a.N:run(name,N,a.dps,a.cycles,a.out)
if __name__=='__main__':main()

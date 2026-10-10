"""Revision-stage synthetic analyses. See PROTOCOL.md; no patient data."""
from pathlib import Path
import sys, json, csv, math, platform, argparse, hashlib
import numpy as np
from scipy.stats import beta
import scipy
from collections import defaultdict
R=Path(__file__).resolve().parent
OLD=R.parents[1]/'extensions/bounded_correction'
EPS=1e-10
PRIMS=('source_gold_plus','source_gold_minus','source_se','source_sp','target_proxy_plus','target_proxy_minus','target_audit_recorded')

def f(r,l,s,b):
    k=l-1
    return r*(1+k*b)/(1+k*((1-r)*s+r*b))

def box(method,recording):
    u,v=(1.,1.) if recording=='accurate' else (.8,.95)
    return ((4.,4.),(u,u),(v,v)) if method=='oracle' else ((3.,5.),(.95,1.) if recording=='accurate' else (.75,.85),(.95,1.) if recording=='accurate' else (.90,1.))

def joint(ci,mechanism,grid=64):
    """Outward interval-cell bounds. No point-grid optimization or sharpness claim."""
    s0,s1=max(0,ci[2][0]-.05),min(1,ci[2][1]+.05)
    b0,b1=max(0,1-ci[3][1]-.05),min(1,1-ci[3][0]+.05)
    sl,su=np.linspace(s0,s1,grid+1)[:-1,None],np.linspace(s0,s1,grid+1)[1:,None]
    bl,bu=np.linspace(b0,b1,grid+1)[None,:-1],np.linspace(b0,b1,grid+1)[None,1:]
    sl=np.maximum(sl,bl+.5); bu=np.minimum(bu,su-.5)
    valid=(sl<=su+EPS)&(bl<=bu+EPS)
    # Invalid cells remain masked; safe denominators avoid warnings there.
    denlo=np.maximum(su-bu,1e-12); denhi=np.maximum(sl-bl,1e-12)
    lower=[];upper=[]
    for g in range(2):
        lower.append(np.maximum(max(0,ci[g][0]-.1),(ci[4+g][0]-bu)/denlo))
        upper.append(np.minimum(min(1,ci[g][1]+.1),(ci[4+g][1]-bl)/denhi))
    lam,u,v=mechanism; a=ci[6]
    if a[0]>u[1]+EPS or a[1]<1-v[1]-EPS:return None
    rl=max(0,(a[0]+v[0]-1)/(u[1]+v[0]-1))
    ru=min(1,(a[1]+v[1]-1)/(u[0]+v[1]-1))
    if rl>ru+EPS:return None
    flo=np.minimum(f(rl,lam[0],su,bl),f(rl,lam[1],su,bl))
    fhi=np.maximum(f(ru,lam[0],sl,bu),f(ru,lam[1],sl,bu))
    lower[0]=np.maximum(lower[0],flo);upper[0]=np.minimum(upper[0],fhi)
    valid &= (lower[0]<=upper[0]+EPS)&(lower[1]<=upper[1]+EPS)
    if not valid.any():return None
    lo=(2*lower[0]-upper[1]-.2)/2.4
    hi=(2*upper[0]-lower[1]-.2)/2.4
    return float(lo[valid].min()-EPS),float(hi[valid].max()+EPS)

def classify(res,truth):
    if res is None:return dict(lower='',upper='',width='',covered=0,correct_sign=0,wrong_sign=0,indeterminate=0,infeasible=1)
    lo,hi=res;pos=lo>EPS;neg=hi<-EPS
    return dict(lower=lo,upper=hi,width=hi-lo,covered=int(lo-EPS<=truth<=hi+EPS),correct_sign=int((pos and truth>EPS)or(neg and truth<-EPS)),wrong_sign=int((pos and truth<=EPS)or(neg and truth>=-EPS)),indeterminate=int(not(pos or neg)),infeasible=0)

def mc(x):
    n=len(x);k=sum(x);p=k/n
    return dict(n=n,count=int(k),rate=p,mcse=math.sqrt(p*(1-p)/n),low=0 if k==0 else float(beta.ppf(.025,k,n-k+1)),high=1 if k==n else float(beta.ppf(.975,k+1,n-k)),one_sided_95_upper=1 if k==n else float(beta.ppf(.95,k+1,n-k)))

def mean(x):
    return dict(n=len(x),mean=float(np.mean(x)) if x else None,mcse=float(np.std(x,ddof=1)/np.sqrt(len(x))) if len(x)>1 else None)

def writecsv(name,rows):
    with (R/name).open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def summarize(rows,keys):
    groups=defaultdict(list)
    for r in rows:groups[tuple(r[k] for k in keys)].append(r)
    out=[]
    for key,rs in groups.items():
        d=dict(zip(keys,key))
        for m in ('covered','correct_sign','wrong_sign','indeterminate','infeasible'):d[m]=mc([r[m] for r in rs])
        for m in ('width','width_reduction','regret'):
            if m in rs[0]:d[m]=mean([r[m] for r in rs if r[m]!=''])
        out.append(d)
    return out

def run_joint():
    records=list(csv.DictReader((OLD/'primitive_counts_main.csv').open(encoding='utf-8')))
    legacy={tuple(r[k] for k in ('seed','cell','rep','method')):r for r in csv.DictReader((OLD/'method_results_main.csv').open(encoding='utf-8'))}
    output=[]; sensitivity=[]; checked=0
    for i,r in enumerate(records):
        ci=[(float(r[p+'_lower']),float(r[p+'_upper'])) for p in PRIMS]
        results={}
        for method in ('oracle','bounded'):
            old=legacy[tuple(r[k] for k in ('seed','cell','rep'))+(method,)]
            outer=joint(ci,box(method,r['recording']))
            if old['lower']=='':outer=None
            elif outer is not None:
                outer=(max(outer[0],float(old['lower'])-EPS),min(outer[1],float(old['upper'])+EPS))
                if outer[0]>outer[1]+EPS:outer=None
            truth=float(r['true_delta']);metrics=classify(outer,truth)
            if r['all_primitives_cover'] in ('True','1'):
                assert metrics['covered']==1,(i,method,outer,truth)
                checked+=1
            results[method]=outer
            output.append({**{k:r[k] for k in ('seed','cell','world','recording','rep')},'method':method,'true_delta':truth,**metrics,'legacy_width':old['width'],'width_reduction':float(old['width'])-metrics['width'] if old['width']!='' and metrics['width']!='' else ''})
            if int(r['rep'])<20:
                grids=[joint(ci,box(method,r['recording']),g) for g in (32,64,128)]
                for a,b in zip(grids,grids[1:]):
                    if b is not None:assert a is not None and a[0]<=b[0]+1e-8 and a[1]>=b[1]-1e-8
                sensitivity.append({**{k:r[k] for k in ('seed','cell','rep')},'method':method,**{f'{g}_{end}':x[j] if x else '' for g,x in zip((32,64,128),grids) for j,end in enumerate(('lower','upper'))}})
        if results['oracle'] and results['bounded']:
            assert results['bounded'][0]<=results['oracle'][0]+1e-8 and results['bounded'][1]>=results['oracle'][1]-1e-8
        if i%1000==0:print('E3 records',i,flush=True)
    writecsv('joint_results.csv',output);writecsv('grid_sensitivity.csv',sensitivity)
    pop=[]
    for w,p,pm,s,sp in [('B',19/60,7/15,.72,.88),('C',.35,.5,.7,.9)]:
        for rec,u,v in [('accurate',1,1),('fallible',.8,.95)]:
            a=(u+v-1)*p*(1+3*s)/(1+3*.31)+1-v
            ci=[(x,x) for x in (.4,.55,.67,.93,.31,.4,a)]
            for m in ('oracle','bounded'):
                for g in (64,128,256):
                    pop.append(dict(world=w,recording=rec,method=m,grid=g,**classify(joint(ci,box(m,rec),g),(2*p-pm-.2)/2.4)))
    writecsv('joint_population.csv',pop)
    summary=dict(evidence_records=len(records),method_records=len(output),primitive_event_containment_checks=checked,summary=summarize(output,['world','recording','method']),per_seed=summarize(output,['seed','world','recording','method']),population=pop)
    (R/'joint_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')

def run_budget():
    output=[];seeds=(2026101001,2026101002,2026101003)
    worlds={'A':(29/60,19/30),'B':(19/60,7/15),'C':(.35,.5)}
    for seed in seeds:
        rng=np.random.default_rng(seed)
        for world,(p,pm) in worlds.items():
            truth=0.0 if world=='C' else (2*p-pm-.2)/2.4
            for n in (200,800,3200,12800):
                for allocation in (1/3,1/2,2/3):
                    np_=int(n*allocation);nm=n-np_
                    kp=rng.binomial(np_,p,1000);km=rng.binomial(nm,pm,1000)
                    def cp(k,n):
                        return np.where(k==0,0,beta.ppf(.0125,k,n-k+1)),np.where(k==n,1,beta.ppf(.9875,k+1,n-k))
                    lp,up=cp(kp,np_);lm,um=cp(km,nm)
                    lo=(2*lp-um-.2)/2.4;hi=(2*up-lm-.2)/2.4
                    for rep in range(1000):
                        choice=1 if max(0,-lo[rep])<=max(0,hi[rep]) else 2
                        regret=max(0,-truth) if choice==1 else max(0,truth)
                        output.append(dict(seed=seed,world=world,budget=n,allocation=allocation,rep=rep,n_plus=np_,n_minus=nm,k_plus=int(kp[rep]),k_minus=int(km[rep]),true_delta=truth,**classify((float(lo[rep]),float(hi[rep])),truth),choice=choice,regret=regret))
    writecsv('budget_results.csv',output)
    (R/'budget_summary.json').write_text(json.dumps(dict(records=len(output),seeds=seeds,summary=summarize(output,['world','budget','allocation']),per_seed=summarize(output,['seed','world','budget','allocation'])),indent=2),encoding='utf-8')
    print('E4 complete',len(output),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--part',choices=['joint','budget','all'],default='all');args=parser.parse_args()
    (R/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,protocol_sha256=hashlib.sha256((R/'PROTOCOL.md').read_bytes()).hexdigest()),indent=2),encoding='utf-8')
    if args.part in ('joint','all'):run_joint()
    if args.part in ('budget','all'):run_budget()

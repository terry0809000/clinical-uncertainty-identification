"""Separate all-row reconstruction: imports no simulation/correction implementation.

Same-family algorithmic audit, not external independent authorship. Uses a separate
reference-corner inversion and reordered LP, solved with highs-ds.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.stats import beta

HERE=Path(__file__).resolve().parent
NAMES=('source_gold_plus','source_gold_minus','source_se','source_sp',
       'target_proxy_plus','target_proxy_minus','target_audit_recorded')
METHODS=('naive','oracle','bounded','wrong_selection_bound')
SEEDS=(2026093001,2026093002,2026093003)
METRICS=('true_covered','correct_sign','wrong_sign','indeterminate','infeasible')
TRUTH={'B':(19/60,7/15,.72,.88,-1/72),'C':(.35,.5,.7,.9,0.)}


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()


@lru_cache(None)
def exact_binomial(x,n):
    if n==0: return (0.,1.)
    return (float(beta.ppf(.05/14,x,n-x+1)) if x else 0.,
            float(beta.ppf(1-.05/14,x+1,n-x)) if x<n else 1.)


def box(method,recording):
    e,f=(1.,1.) if recording=='accurate' else (.8,.95)
    if method=='naive': return (1.,1.),(1.,1.),(1.,1.)
    if method=='oracle': return (4.,4.),(e,e),(f,f)
    if method=='wrong_selection_bound': return (1.,1.25),(e,e),(f,f)
    return (3.,5.), ((.95,1.) if recording=='accurate' else (.75,.85)), \
        ((.95,1.) if recording=='accurate' else (.9,1.))


def inverse_by_corners(a,q,se,b):
    lam,us,vs=b
    # All reference-box corners plus intersections with r=0 and r=1 boundaries.
    uc=sorted({x for x in (*us,*a) if us[0]<=x<=us[1]})
    vc=sorted({x for x in (*vs,1-a[0],1-a[1]) if vs[0]<=x<=vs[1]})
    risks=[]
    for u,v in itertools.product(uc,vc):
        al,ah=max(a[0],1-v),min(a[1],u)
        if al<=ah:
            risks.extend(((al+v-1)/(u+v-1),(ah+v-1)/(u+v-1)))
    if not risks: return None
    factors=[(1+(l-1)*qq)/(1+(l-1)*ss) for l,qq,ss in itertools.product(lam,q,se)]
    lo,hi=max(0.,min(risks)*min(factors)),min(1.,max(risks)*max(factors))
    return None if lo>hi else (lo,hi)


def independent_lp(cis,p):
    if p is None: return None
    # Variable order p_plus, p_minus, reciprocal Youden, false-positive/Youden.
    lower_s,upper_s=max(0.,cis[2][0]-.05),min(1.,cis[2][1]+.05)
    lower_c,upper_c=max(0.,cis[3][0]-.05),min(1.,cis[3][1]+.05)
    pl=[max(0.,cis[g][0]-.1) for g in (0,1)]
    pu=[min(1.,cis[g][1]+.1) for g in (0,1)]
    pl[0],pu[0]=max(pl[0],p[0]),min(pu[0],p[1])
    if any(x>y for x,y in zip(pl,pu)): return None
    aa=[];bb=[]
    # qL z <= p_g+w <= qU z
    for g in (0,1):
        a=[0.,0.,-cis[4+g][1],1.];a[g]=1.
        aa.append(a);bb.append(0.)
        a=[0.,0.,cis[4+g][0],-1.];a[g]=-1.
        aa.append(a);bb.append(0.)
    # Se=(1+w)/z and Sp=1-w/z.
    aa.extend(([0,0,-upper_s,1],[0,0,lower_s,-1],
               [0,0,-(1-lower_c),1],[0,0,1-upper_c,-1]))
    bb.extend((-1,1,0,0))
    obj=np.array([5/6,-5/12,0.,0.])
    bounds=((pl[0],pu[0]),(pl[1],pu[1]),(1.,2.),(None,None))
    rlo=linprog(obj,A_ub=aa,b_ub=bb,bounds=bounds,method='highs-ds')
    if rlo.status==2: return None
    assert rlo.success,rlo.message
    rhi=linprog(-obj,A_ub=aa,b_ub=bb,bounds=bounds,method='highs-ds')
    assert rhi.success,rhi.message
    return float(rlo.fun-1/12),float(-rhi.fun-1/12)


def flags(interval,truth):
    if interval is None: return dict(zip(MR:=(METRICS), (0,0,0,0,1)))
    lo,hi=interval
    pos,neg=lo>1e-10,hi< -1e-10
    return {'true_covered':int(lo-1e-10<=truth<=hi+1e-10),
        'correct_sign':int((pos and truth>0) or (neg and truth<0)),
        'wrong_sign':int((pos and truth<=0) or (neg and truth>=0)),
        'indeterminate':int(not(pos or neg)),'infeasible':0}


def main():
    summary=json.loads((HERE/'summary_main.json').read_text())
    evidence=list(csv.DictReader((HERE/'primitive_counts_main.csv').open(newline='',encoding='utf-8')))
    results=list(csv.DictReader((HERE/'method_results_main.csv').open(newline='',encoding='utf-8')))
    assert len(evidence)==12000 and len(results)==48000
    emap={}
    maxdiff=defaultdict(float)
    def close(label,x,y,tol=1e-10):
        d=abs(float(x)-float(y));maxdiff[label]=max(maxdiff[label],d)
        assert d<=tol,(label,x,y,d)
    for row in evidence:
        seed,cell,rep=int(row['seed']),int(row['cell']),int(row['rep'])
        key=(seed,cell,rep)
        assert key not in emap
        world='B' if cell<2 else 'C'
        rec='accurate' if cell%2==0 else 'fallible'
        assert row['world']==world and row['recording']==rec
        p,pm,s,sp,d=TRUTH[world]
        u,v=(1.,1.) if rec=='accurate' else (.8,.95)
        sel=p*(1+3*s)/(1+3*.31);a=(u+v-1)*sel+1-v
        for name,value in {'true_delta':d,'lambda_true':4,'ref_se_true':u,'ref_sp_true':v,
            'p_plus_true':p,'p_minus_true':pm,'target_se_true':s,'target_sp_true':sp,
            'selected_true_risk':sel,'audit_recorded_probability':a}.items():
            close('truth_metadata',row[name],value)
        cis=[];ks=[];ns=[]
        for name in NAMES:
            k,n=int(row[name+'_success']),int(row[name+'_denominator'])
            assert 0<=k<=n
            ci=exact_binomial(k,n);cis.append(ci);ks.append(k);ns.append(n)
            for suffix,val in zip(('lower','upper'),ci):close('cp',row[name+'_'+suffix],val,1e-12)
        assert ns==(expected:=[1000,1000,ks[0]+ks[1],2000-ks[0]-ks[1],1000,1000,800]),(ns,expected)
        covered=all(lo<=t<=hi for (lo,hi),t in zip(cis,(.4,.55,.67,.93,.31,.4,a)))
        assert int(row['all_primitives_cover'])==covered
        emap[key]=(cis,world,rec,d,p,covered)
    expected=set(itertools.product(SEEDS,range(4),range(1000)))
    assert set(emap)==expected
    counts=defaultdict(lambda:defaultdict(int));widths=defaultdict(list)
    methods_by_key=defaultdict(dict);seen=set();valid_events=0
    for i,row in enumerate(results):
        seed,cell,rep=int(row['seed']),int(row['cell']),int(row['rep'])
        method=row['method'];key=(seed,cell,rep)
        assert method in METHODS and (*key,method) not in seen
        seen.add((*key,method))
        cis,world,rec,truth,p,covered=emap[key]
        assert row['world']==world and row['recording']==rec
        s=(max(0.,cis[2][0]-.05),min(1.,cis[2][1]+.05))
        pb=inverse_by_corners(cis[6],cis[4],s,box(method,rec))
        if pb is None: assert row['p_lower']==row['p_upper']==''
        else:
            close('inverse',pb[0],row['p_lower']);close('inverse',pb[1],row['p_upper'])
        interval=independent_lp(cis,pb)
        expected_flags=flags(interval,truth)
        for met in METRICS: assert int(row[met])==expected_flags[met],(key,method,met)
        if interval is None: assert row['lower']==row['upper']==row['width']==''
        else:
            close('lp',interval[0],row['lower']);close('lp',interval[1],row['upper'])
            close('width',interval[1]-interval[0],row['width'])
        if covered and method in ('oracle','bounded'):
            assert pb is not None and pb[0]-1e-10<=p<=pb[1]+1e-10
            assert expected_flags['true_covered'] and not expected_flags['wrong_sign']
            valid_events+=1
        methods_by_key[key][method]=(pb,interval)
        for group in ((seed,world,rec,method),('pooled',world,rec,method)):
            counts[group]['n']+=1
            for met in METRICS:counts[group][met]+=expected_flags[met]
            if interval is not None:widths[group].append(interval[1]-interval[0])
        if (i+1)%8000==0:print(f'Audited {i+1}/48000 method rows',flush=True)
    assert len(seen)==48000
    pairs=defaultdict(list)
    for key,methods in methods_by_key.items():
        assert set(methods)==set(METHODS)
        ob,oe=methods['oracle'];bb,be=methods['bounded']
        if ob is not None: assert bb is not None and bb[0]<=ob[0]+1e-10 and bb[1]>=ob[1]-1e-10
        if oe is not None:
            assert be is not None and be[0]<=oe[0]+1e-9 and be[1]>=oe[1]-1e-9
            diff=(be[1]-be[0])-(oe[1]-oe[0])
            seed,cell,_=key;_,world,rec,_,_,_=emap[key]
            for group in ((seed,world,rec),('pooled',world,rec)):pairs[group].append(diff)
    for cell in summary['cells']:
        key=(cell['seed'],cell['world'],cell['recording'],cell['method']);c=counts[key];n=c['n']
        for met in METRICS:
            k=c[met];m=cell['metrics'][met];assert m['count']==k and m['n']==n
            close('summaries',m['proportion'],k/n)
            close('summaries',m['mcse'],math.sqrt((k/n)*(1-k/n)/n))
            bounds=(0. if k==0 else beta.ppf(.025,k,n-k+1),1. if k==n else beta.ppf(.975,k+1,n-k))
            for actual,expected in zip(m['mc_95_interval'],bounds):close('mc_interval',actual,expected)
            close('mc_interval',m['one_sided_95_upper'],1. if k==n else beta.ppf(.95,k+1,n-k))
        assert sum(c[m] for m in ('correct_sign','wrong_sign','indeterminate','infeasible'))==n
        if widths[key]:close('mean_width',cell['mean_width_feasible_only'],np.mean(widths[key]))
        else:assert cell['mean_width_feasible_only'] is None
    for cell in summary['paired_width_bounded_minus_oracle']:
        data=pairs[(cell['seed'],cell['world'],cell['recording'])]
        assert cell['n_both_feasible']==len(data)
        for name,value in {'mean':np.mean(data),'mcse':np.std(data,ddof=1)/math.sqrt(len(data)),
                           'minimum':min(data),'maximum':max(data)}.items():close('paired_width',cell[name],value)
    assert valid_events==summary['valid_primitive_event_containment_checks']
    empirical_failures=[]
    for (seed,world,rec,method),c in counts.items():
        if seed=='pooled' and method in ('oracle','bounded'):
            k,n=c['true_covered'],c['n']
            upper=1. if k==n else beta.ppf(1-.05/8,k+1,n-k)
            if upper<.95:empirical_failures.append((world,rec,method,k,n,float(upper)))
    assert not empirical_failures
    for name,value in summary['sha256'].items():
        path=HERE/name
        if not path.exists():path=HERE.parents[2]/'review_packages/resubmission_20260930/companion_original/simulation'/name
        assert sha(path)==value,(name,'hash mismatch')
    result={'status':'PASS_DETERMINISTIC','review_independence':'same_family_provisional; separate implementation, no imports from analysis',
        'n_evidence_records':12000,'n_method_rows':48000,'seeds':SEEDS,'n_per_seed_generating_cell':1000,
        'max_abs_discrepancy':dict(maxdiff),'valid_primitive_event_checks':valid_events,
        'bounded_contains_oracle_rows':12000,'empirical_valid_coverage_failures':empirical_failures,
        'verified':['complete unique keys','all counts/denominators and CP intervals','generator metadata',
                    'corner-enumerated reference inverse','all independently reordered LPs via highs-ds',
                    'all row flags, summaries, MC intervals and paired widths','all-row oracle nesting',
                    'simultaneous primitive-event containment','input/script/protocol hashes'],
        'sha256':{p.name:sha(p) for p in (HERE/'PROTOCOL.md',HERE/'run_bounded_correction.py',Path(__file__),
                                      HERE/'primitive_counts_main.csv',HERE/'method_results_main.csv',HERE/'summary_main.json')}}
    (HERE/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('PASS_DETERMINISTIC',json.dumps(dict(maxdiff)),flush=True)


if __name__=='__main__':main()

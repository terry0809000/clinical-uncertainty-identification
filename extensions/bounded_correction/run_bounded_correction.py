"""Frozen synthetic bounded audit correction; see PROTOCOL.md. CPU only."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import platform
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1] # Portable companion root; scientific algorithm unchanged.
ORIGINAL = PROJECT / 'legacy/simulation/run_study.py'
# Reuse the archived LP without editing it; the independent checker does not import it.
import sys
spec = importlib.util.spec_from_file_location('archived_study', ORIGINAL)
old = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = old
spec.loader.exec_module(old)

SEEDS = (2026093001, 2026093002, 2026093003)
PRIMITIVES = ('source_gold_plus', 'source_gold_minus', 'source_se', 'source_sp',
              'target_proxy_plus', 'target_proxy_minus', 'target_audit_recorded')
METHODS = ('naive', 'oracle', 'bounded', 'wrong_selection_bound')
METRICS = ('true_covered', 'correct_sign', 'wrong_sign', 'indeterminate', 'infeasible')
WORLDS = {'B': (19 / 60, 7 / 15, .72, .88, -1 / 72),
          'C': (.35, .50, .70, .90, 0.0)}
DESIGNS = tuple((w, rec) for w in WORLDS for rec in ('accurate', 'fallible'))
PROTOCOL_SHA = '9E734121754BDAECCF2B67DE18A8E71F9FD64530BF59C39EE52B8A40B39AAFFD'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()


@lru_cache(maxsize=None)
def cp(k, n):
    if n == 0:
        return (0.0, 1.0)
    tail = .05 / 14
    return (0.0 if k == 0 else float(beta.ppf(tail, k, n-k+1)),
            1.0 if k == n else float(beta.ppf(1-tail, k+1, n-k)))


def mechanism(method, recording):
    u, v = (1., 1.) if recording == 'accurate' else (.80, .95)
    if method == 'naive':
        return ((1., 1.), (1., 1.), (1., 1.))
    if method == 'oracle':
        return ((4., 4.), (u, u), (v, v))
    if method == 'wrong_selection_bound':
        return ((1., 1.25), (u, u), (v, v))
    return ((3., 5.), (.95, 1.) if recording == 'accurate' else (.75, .85),
            (.95, 1.) if recording == 'accurate' else (.90, 1.))


def inverse(a, q, se, box):
    lam, u, v = box
    assert lam[0] >= 1 and u[0] + v[0] > 1
    if a[0] > u[1] or a[1] < 1-v[1]:
        return None
    rl = max(0., (a[0]+v[0]-1)/(u[1]+v[0]-1))
    ru = min(1., (a[1]+v[1]-1)/(u[0]+v[1]-1))
    factors_l = [(1+(l-1)*q[0])/(1+(l-1)*se[1]) for l in lam]
    factors_u = [(1+(l-1)*q[1])/(1+(l-1)*se[0]) for l in lam]
    lo, hi = max(0., rl*min(factors_l)), min(1., ru*max(factors_u))
    return None if lo > hi else (lo, hi)


def projection(cis, bound):
    if bound is None:
        return None
    return old.envelope(cis[:2], cis[2], cis[3], cis[4:6], bound)


def classify(result, truth):
    if result is None:
        return {'lower': '', 'upper': '', 'width': '', 'true_covered': 0,
                'correct_sign': 0, 'wrong_sign': 0, 'indeterminate': 0, 'infeasible': 1}
    lo, hi = result
    positive, negative = lo > 1e-10, hi < -1e-10
    wrong = (positive and truth <= 0) or (negative and truth >= 0)
    correct = (positive and truth > 0) or (negative and truth < 0)
    return {'lower': lo, 'upper': hi, 'width': hi-lo,
            'true_covered': int(lo-1e-10 <= truth <= hi+1e-10),
            'correct_sign': int(correct), 'wrong_sign': int(wrong),
            'indeterminate': int(not (positive or negative)), 'infeasible': 0}


def deterministic_checks():
    cases = 0
    # Verify against forward laws, including reversed q versus Se ordering.
    for p in (0., .1, .35, .8, 1.):
        for s in (.1, .7, 1.):
            for specificity in (.1, .9, 1.):
                q = s*p+(1-specificity)*(1-p)
                for lam in (1., 3., 5.):
                    r = p*(1+(lam-1)*s)/(1+(lam-1)*q)
                    for u,v in ((1.,1.),(.8,.95),(.6,.6)):
                        a = (u+v-1)*r+1-v
                        ab = (max(0., a-.01), min(1., a+.01))
                        qb = (max(0., q-.01), min(1., q+.01))
                        sb = (max(0., s-.01), min(1., s+.01))
                        ub = (max(.51,u-.03), min(1.,u+.03))
                        vb = (max(.51,v-.03), min(1.,v+.03))
                        result = inverse(ab,qb,sb,((max(1.,lam-.5),lam+.5),ub,vb))
                        assert result is not None and result[0]-1e-12 <= p <= result[1]+1e-12
                        identity = inverse(ab,qb,sb,((1.,1.),(1.,1.),(1.,1.)))
                        assert np.allclose(identity, ab, atol=1e-12, rtol=0)
                        cases += 1
    assert inverse((.91,.99),(.3,.4),(.6,.8),((3.,5.),(.75,.85),(.9,1.))) is None
    assert inverse((0.,.01),(.3,.4),(.6,.8),((3.,5.),(.75,.85),(.8,.9))) is None
    return {'forward_containment_and_identity_cases': cases, 'physical_impossibility_cases': 2}


def population_checks():
    rows = []
    for world, recording in DESIGNS:
        p, pm, s, sp, truth = WORLDS[world]
        u,v = (1.,1.) if recording == 'accurate' else (.8,.95)
        q = s*p+(1-sp)*(1-p)
        assert abs(q-.31) < 1e-12 and abs(s*pm+(1-sp)*(1-pm)-.4)<1e-12
        r = p*(1+3*s)/(1+3*q)
        a = (u+v-1)*r+1-v
        cis = tuple((x,x) for x in (.4,.55,.67,.93,.31,.4,a))
        for method in METHODS:
            bound = inverse(cis[6],cis[4],(.62,.72),mechanism(method,recording))
            out = projection(cis,bound)
            rows.append({'world':world,'recording':recording,'method':method,
                         'a':a,'selected_risk':r,'inverse':bound,
                         'outer_envelope':out,'sharp':False,**classify(out,truth)})
    return rows


def mc(k,n):
    p = k/n
    return {'count': k, 'n': n, 'proportion': p, 'mcse': math.sqrt(p*(1-p)/n),
            'mc_95_interval': [0. if k==0 else float(beta.ppf(.025,k,n-k+1)),
                               1. if k==n else float(beta.ppf(.975,k+1,n-k))],
            'one_sided_95_upper': 1. if k==n else float(beta.ppf(.95,k+1,n-k))}


def main(reps, label):
    assert sha(HERE/'PROTOCOL.md') == PROTOCOL_SHA
    checks = deterministic_checks()
    pop = population_checks()
    counts = defaultdict(lambda:defaultdict(int))
    widths = defaultdict(list)
    paired_width = defaultdict(list)
    valid_invariants = 0
    raw = HERE/f'primitive_counts_{label}.csv'
    outputs = HERE/f'method_results_{label}.csv'
    summary_path = HERE/f'summary_{label}.json'
    base = ['seed','cell','world','recording','rep','true_delta','lambda_true','ref_se_true',
            'ref_sp_true','p_plus_true','p_minus_true','target_se_true','target_sp_true',
            'selected_true_risk','audit_recorded_probability','all_primitives_cover']
    raw_fields = base + [f'{p}_{x}' for p in PRIMITIVES for x in ('success','denominator','lower','upper')]
    out_fields = ['seed','cell','world','recording','rep','method','p_lower','p_upper',
                  'lower','upper','width']+list(METRICS)
    with raw.open('w',newline='',encoding='utf-8') as fraw, outputs.open('w',newline='',encoding='utf-8') as fout:
        rw,ow = csv.DictWriter(fraw,fieldnames=raw_fields),csv.DictWriter(fout,fieldnames=out_fields)
        rw.writeheader(); ow.writeheader()
        for seed in SEEDS:
            children=np.random.SeedSequence(seed).spawn(len(DESIGNS))
            for cell,(world,recording) in enumerate(DESIGNS):
                rng=np.random.default_rng(children[cell])
                p,pm,s,sp,truth=WORLDS[world]
                u,v=(1.,1.) if recording=='accurate' else (.8,.95)
                selected=p*(1+3*s)/(1+3*.31)
                a=(u+v-1)*selected+1-v
                truths=(.4,.55,.67,.93,.31,.4,a)
                for rep in range(reps):
                    k=[int(rng.binomial(1000,pr)) for pr in (.4,.55)]
                    sy=sum(int(rng.binomial(kg,.67)) for kg in k)
                    sn=sum(int(rng.binomial(1000-kg,.93)) for kg in k)
                    qp,qm=int(rng.binomial(1000,.31)),int(rng.binomial(1000,.4))
                    ar=int(rng.binomial(800,a))
                    obs=((k[0],1000),(k[1],1000),(sy,sum(k)),(sn,2000-sum(k)),(qp,1000),(qm,1000),(ar,800))
                    cis=tuple(cp(x,n) for x,n in obs)
                    all_cover=all(lo<=t<=hi for (lo,hi),t in zip(cis,truths))
                    keybase={'seed':seed,'cell':cell,'world':world,'recording':recording,'rep':rep}
                    rw.writerow({**keybase,'true_delta':truth,'lambda_true':4.,'ref_se_true':u,'ref_sp_true':v,
                        'p_plus_true':p,'p_minus_true':pm,'target_se_true':s,'target_sp_true':sp,
                        'selected_true_risk':selected,'audit_recorded_probability':a,'all_primitives_cover':int(all_cover),
                        **{f'{name}_{part}':value for name,(x,n),(lo,hi) in zip(PRIMITIVES,obs,cis)
                           for part,value in (('success',x),('denominator',n),('lower',lo),('upper',hi))}})
                    target_se=(max(0.,cis[2][0]-.05),min(1.,cis[2][1]+.05))
                    record={}
                    for method in METHODS:
                        bound=inverse(cis[6],cis[4],target_se,mechanism(method,recording))
                        interval=projection(cis,bound)
                        result=classify(interval,truth)
                        record[method]=(bound,interval,result)
                        if method in ('oracle','bounded') and all_cover:
                            assert bound is not None and bound[0]-1e-10<=p<=bound[1]+1e-10
                            assert result['true_covered'] and not result['wrong_sign']
                            valid_invariants+=1
                        ow.writerow({**keybase,'method':method,'p_lower':'' if bound is None else bound[0],
                                     'p_upper':'' if bound is None else bound[1],**result})
                        for group in ((seed,world,recording,method),('pooled',world,recording,method)):
                            counts[group]['n']+=1
                            for met in METRICS: counts[group][met]+=result[met]
                            if not result['infeasible']: widths[group].append(result['width'])
                    ib,io=record['bounded'][0],record['oracle'][0]
                    if io is not None:
                        assert ib is not None and ib[0]<=io[0]+1e-10 and ib[1]>=io[1]-1e-10
                    eb,eo=record['bounded'][1],record['oracle'][1]
                    if eo is not None:
                        assert eb is not None and eb[0]<=eo[0]+1e-9 and eb[1]>=eo[1]-1e-9
                        difference=(eb[1]-eb[0])-(eo[1]-eo[0])
                        for group in ((seed,world,recording),('pooled',world,recording)):
                            paired_width[group].append(difference)
                print(f'{label}: seed {seed}, cell {cell}, {world}/{recording}, {reps} records complete',flush=True)
    cells=[]
    for key,c in counts.items():
        n=c['n']
        assert sum(c[m] for m in ('correct_sign','wrong_sign','indeterminate','infeasible'))==n
        cells.append({'seed':key[0],'world':key[1],'recording':key[2],'method':key[3],
                      'metrics':{met:mc(c[met],n) for met in METRICS},
                      'mean_width_feasible_only':float(np.mean(widths[key])) if widths[key] else None})
    summary={'status':'COMPLETE_PENDING_INDEPENDENT_AUDIT','same_family':'provisional',
        'label':label,'repetitions_per_seed_cell':reps,'seeds':SEEDS,'n_evidence_records':len(SEEDS)*4*reps,
        'n_method_rows':len(SEEDS)*4*reps*4,'protocol_sha256':PROTOCOL_SHA,
        'deterministic_checks':checks,'valid_primitive_event_containment_checks':valid_invariants,
        'population_outer_envelopes_not_sharp':pop,'cells':cells,
        'paired_width_bounded_minus_oracle':[{'seed':k[0],'world':k[1],'recording':k[2],
            'n_both_feasible':len(v),'mean':float(np.mean(v)),
            'mcse':float(np.std(v,ddof=1)/math.sqrt(len(v))) if len(v)>1 else None,
            'minimum':min(v),'maximum':max(v)} for k,v in paired_width.items()],
        'sha256':{p.name:sha(p) for p in (HERE/'PROTOCOL.md',Path(__file__),ORIGINAL,raw,outputs)},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'limits':['synthetic externally stipulated mechanism ranges','common target proxy channel retained',
                  'constant reference channel conditional on true outcome','outer projection not sharp set',
                  'coverage simulation only of true point','no efficiency or clinical validity claim']}
    summary_path.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(summary_path,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--reps',type=int,choices=(10,1000),required=True)
    parser.add_argument('--label',choices=('pilot','main'),required=True)
    args=parser.parse_args()
    assert (args.reps,args.label) in ((10,'pilot'),(1000,'main'))
    main(args.reps,args.label)

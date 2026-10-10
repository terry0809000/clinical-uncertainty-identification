"""Separate reconstruction of E3/E4; imports no author analysis functions."""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import csv
import json
import math
import time
import numpy as np
from scipy.special import betaincinv

ROOT = Path(__file__).resolve().parents[1]
NEW = ROOT / 'experiments'
OLD = ROOT.parent / 'extensions/bounded_correction'
EPS = 1e-10
started = time.time()
max_error = defaultdict(float)
checks = defaultdict(int)

def rows(path):
    with path.open(encoding='utf-8') as handle:
        return list(csv.DictReader(handle))

def close(label, observed, expected, tol=2e-9):
    discrepancy = float(np.max(np.abs(np.asarray(observed, float)-np.asarray(expected, float))))
    max_error[label] = max(max_error[label], discrepancy)
    assert discrepancy <= tol, (label, discrepancy, observed, expected)

def cp(success, total, error):
    k = np.asarray(success, float)
    n = np.asarray(total, float)
    lo = np.zeros_like(k)
    hi = np.ones_like(k)
    left = (k > 0) & (n > 0)
    right = (k < n) & (n > 0)
    lo[left] = betaincinv(k[left], n[left]-k[left]+1, error/2)
    hi[right] = betaincinv(k[right]+1, n[right]-k[right], 1-error/2)
    return np.stack((lo, hi), axis=-1)

def ref_interval(a, u, v):
    # Enumerate the physical polytope's boundary candidates rather than import inversion.
    possibilities = []
    for aa in a:
        for uu in u:
            for vv in v:
                rr = (aa+vv-1)/(uu+vv-1)
                if -1e-12 <= rr <= 1+1e-12:
                    possibilities.append(min(1., max(0., rr)))
    if max(a[0], 1-v[1]) <= min(a[1], 1-v[0]):
        possibilities.append(0.)
    if max(a[0], u[0]) <= min(a[1], u[1]):
        possibilities.append(1.)
    if not possibilities:
        return None
    return min(possibilities), max(possibilities)

def mechanism(method, recording):
    if method == 'oracle':
        return (4., 4.), ((1.,1.) if recording=='accurate' else (.8,.8)), ((1.,1.) if recording=='accurate' else (.95,.95))
    return (3.,5.), ((.95,1.) if recording=='accurate' else (.75,.85)), ((.95,1.) if recording=='accurate' else (.9,1.))

def inverse_odds(r, relative_weight, sensitivity, false_positive):
    # Selection multiplies disease odds by (1+k*s)/(1+k*b).
    diseased_weight = 1+(relative_weight-1)*sensitivity
    nondiseased_weight = 1+(relative_weight-1)*false_positive
    return r*nondiseased_weight / ((1-r)*diseased_weight+r*nondiseased_weight)

def reconstruct_joint(intervals, spec, count=64):
    lam, u, v = spec
    r = ref_interval(intervals[6], u, v)
    if r is None:
        return None
    s_limits = max(0.,intervals[2,0]-.05), min(1.,intervals[2,1]+.05)
    b_limits = max(0.,1-intervals[3,1]-.05), min(1.,1-intervals[3,0]+.05)
    s_grid = np.linspace(*s_limits, count+1)
    b_grid = np.linspace(*b_limits, count+1)
    # Flatten in the opposite axis order to the analysis implementation.
    smin = np.tile(s_grid[:-1], count)
    smax = np.tile(s_grid[1:], count)
    bmin = np.repeat(b_grid[:-1], count)
    bmax = np.repeat(b_grid[1:], count)
    keep = smax-bmin >= .5-EPS
    smin = np.maximum(smin, bmin+.5)
    bmax = np.minimum(bmax, smax-.5)
    keep &= (smin <= smax+EPS) & (bmin <= bmax+EPS)
    ilower = np.empty((2,count*count))
    iupper = np.empty_like(ilower)
    for g in (0,1):
        # Invert q <= bmax+(smax-bmax)*p and q >= bmin+(smin-bmin)*p.
        l = (intervals[g+4,0]-bmax)/np.maximum(smax-bmax,1e-12)
        h = (intervals[g+4,1]-bmin)/np.maximum(smin-bmin,1e-12)
        ilower[g] = np.maximum(l, max(0.,intervals[g,0]-.1))
        iupper[g] = np.minimum(h, min(1.,intervals[g,1]+.1))
    pl = np.minimum.reduce([inverse_odds(r[0], z, smax, bmin) for z in lam])
    ph = np.maximum.reduce([inverse_odds(r[1], z, smin, bmax) for z in lam])
    ilower[0] = np.maximum(ilower[0],pl)
    iupper[0] = np.minimum(iupper[0],ph)
    keep &= np.all(ilower <= iupper+EPS,axis=0)
    if not np.any(keep):
        return None
    lo = (10*ilower[0,keep]-5*iupper[1,keep]-1)/12
    hi = (10*iupper[0,keep]-5*ilower[1,keep]-1)/12
    return float(min(lo)-EPS), float(max(hi)+EPS)

def expected_flags(bounds, truth):
    if bounds is None:
        return dict(covered=0,correct_sign=0,wrong_sign=0,indeterminate=0,infeasible=1)
    lo, hi = bounds
    sign = 1 if lo > EPS else -1 if hi < -EPS else 0
    true_sign = 1 if truth > EPS else -1 if truth < -EPS else 0
    return dict(covered=int(lo-EPS<=truth<=hi+EPS), correct_sign=int(sign!=0 and sign==true_sign),
                wrong_sign=int(sign!=0 and sign!=true_sign), indeterminate=int(sign==0), infeasible=0)

def verify_flags(row,bounds,truth):
    flags = expected_flags(bounds,truth)
    for field,value in flags.items():
        assert int(row[field]) == value, (field,row,value)
    checks['row_flags'] += 1

def summary_check(data, summary, keys):
    groups=defaultdict(list)
    def identity(row):
        return tuple(str(row[k]) for k in keys)
    for row in data:
        groups[identity(row)].append(row)
    assert len(summary)==len(groups)
    for entry in summary:
        rs=groups[identity(entry)]
        for field in ('covered','correct_sign','wrong_sign','indeterminate','infeasible'):
            k=sum(int(r[field]) for r in rs); n=len(rs); p=k/n
            lo,hi=cp(np.array([k]),np.array([n]),.05)[0]
            upper=1. if k==n else float(betaincinv(k+1,n-k,.95))
            vals=dict(n=n,count=k,rate=p,mcse=math.sqrt(p*(1-p)/n),low=lo,high=hi,one_sided_95_upper=upper)
            for metric,value in vals.items():
                close('summary_rate',entry[field][metric],value,2e-12)
        for field in ('width','width_reduction','regret'):
            if field in entry:
                vals=np.array([float(r[field]) for r in rs if r[field]!=''])
                close('summary_mean',entry[field]['mean'],np.mean(vals),2e-12)
                close('summary_mcse',entry[field]['mcse'],np.std(vals,ddof=1)/np.sqrt(len(vals)),2e-12)
                assert entry[field]['n']==len(vals)
        checks['summary_cells']+=1

raw = rows(OLD/'primitive_counts_main.csv')
legacy_rows = rows(OLD/'method_results_main.csv')
legacy = {(r['seed'],r['cell'],r['rep'],r['method']):r for r in legacy_rows}
out = rows(NEW/'joint_results.csv')
assert len(raw)==12000 and len(out)==24000
outmap = {(r['seed'],r['cell'],r['rep'],r['method']):r for r in out}
assert len(outmap)==len(out)
prims = ('source_gold_plus','source_gold_minus','source_se','source_sp','target_proxy_plus','target_proxy_minus','target_audit_recorded')
all_ci = np.empty((len(raw),7,2))
for j,name in enumerate(prims):
    k=np.array([int(r[name+'_success']) for r in raw])
    n=np.array([int(r[name+'_denominator']) for r in raw])
    calculated=cp(k,n,.05/7)
    stored=np.array([[float(r[name+'_lower']),float(r[name+'_upper'])] for r in raw])
    close('primitive_intervals',stored,calculated,2e-12)
    all_ci[:,j]=calculated
    checks['primitive_intervals']+=len(raw)
rawmap={}
for i,r in enumerate(raw):
    key=(r['seed'],r['cell'],r['rep']);rawmap[key]=(r,all_ci[i])
    outputs={}
    for method in ('oracle','bounded'):
        result=reconstruct_joint(all_ci[i],mechanism(method,r['recording']))
        previous=legacy[key+(method,)]
        if previous['lower']=='':
            result=None
        elif result is not None:
            result=(max(result[0],float(previous['lower'])-EPS),min(result[1],float(previous['upper'])+EPS))
            if result[0]>result[1]+EPS:result=None
        observed=outmap[key+(method,)]
        assert result is not None and observed['lower']!=''
        close('joint_endpoints',[observed['lower'],observed['upper']],result)
        close('joint_width',observed['width'],result[1]-result[0])
        close('paired_width_reduction',observed['width_reduction'],float(previous['width'])-(result[1]-result[0]))
        truth=float(r['true_delta']);verify_flags(observed,result,truth)
        if int(r['all_primitives_cover']):
            assert result[0]-EPS<=truth<=result[1]+EPS
            checks['primitive_event_true_containment']+=1
        assert result[0]>=float(previous['lower'])-2*EPS and result[1]<=float(previous['upper'])+2*EPS
        outputs[method]=result
        checks['joint_rows']+=1
    assert outputs['bounded'][0]<=outputs['oracle'][0]+2e-9 and outputs['bounded'][1]>=outputs['oracle'][1]-2e-9
    checks['bounded_oracle_nesting']+=1
    if i%2000==0:print('Independent joint rows',2*i,flush=True)

gridrows=rows(NEW/'grid_sensitivity.csv')
assert len(gridrows)==480
for row in gridrows:
    original,ci=rawmap[(row['seed'],row['cell'],row['rep'])]
    last=None
    for g in (32,64,128):
        interval=reconstruct_joint(ci,mechanism(row['method'],original['recording']),g)
        assert interval is not None
        close('grid_reconstruction',[row[f'{g}_lower'],row[f'{g}_upper']],interval)
        if last is not None:
            assert last[0]<=interval[0]+2e-9 and last[1]>=interval[1]-2e-9
        last=interval
        checks['grid_rows']+=1

pop=rows(NEW/'joint_population.csv')
assert len(pop)==24
for row in pop:
    p,pm,s=(19/60,7/15,.72) if row['world']=='B' else (.35,.5,.7)
    u,v=(1.,1.) if row['recording']=='accurate' else (.8,.95)
    aa=(u+v-1)*p*(1+3*s)/(1+3*.31)+1-v
    primitive=np.array([[x,x] for x in (.4,.55,.67,.93,.31,.4,aa)])
    interval=reconstruct_joint(primitive,mechanism(row['method'],row['recording']),int(row['grid']))
    close('population_reconstruction',[row['lower'],row['upper']],interval)
    truth=-1/72 if row['world']=='B' else 0.
    verify_flags(row,interval,truth)
    assert interval[0] <= truth <= interval[1]
    if row['world']=='B' and row['method']=='oracle':assert interval[1]<0
    checks['population_rows']+=1

# Exact algebraic oracle singleton: lower and upper risk inequalities coincide.
rstar=Fraction(1501,2895); K=rstar*(1+3*Fraction(31,100))
lower=K/(1+3*Fraction(18,25))
upper=(K/3+Fraction(3,25)-Fraction(31,100))/(Fraction(1,3)+Fraction(3,25))
assert lower==upper==Fraction(19,60)
assert (10*lower-5*Fraction(7,15)-1)/12==Fraction(-1,72)
checks['exact_oracle_singleton']=1

budget=rows(NEW/'budget_results.csv')
assert len(budget)==108000
unique={(r['seed'],r['world'],r['budget'],r['allocation'],r['rep']) for r in budget}
assert len(unique)==len(budget)
kp=np.array([int(r['k_plus']) for r in budget]); np_=np.array([int(r['n_plus']) for r in budget])
km=np.array([int(r['k_minus']) for r in budget]); nm=np.array([int(r['n_minus']) for r in budget])
ip=cp(kp,np_,.025);im=cp(km,nm,.025)
low=(10*ip[:,0]-5*im[:,1]-1)/12
high=(10*ip[:,1]-5*im[:,0]-1)/12
close('budget_endpoints',[[r['lower'],r['upper']] for r in budget],np.column_stack((low,high)),2e-12)
for i,row in enumerate(budget):
    truth={'A':1/18,'B':-1/72,'C':0.}[row['world']]
    close('budget_truth',row['true_delta'],truth,2e-15)
    verify_flags(row,(low[i],high[i]),truth)
    choice=1 if max(0,-low[i])<=max(0,high[i]) else 2
    loss=max(0,-truth) if choice==1 else max(0,truth)
    assert int(row['choice'])==choice
    close('budget_regret',row['regret'],loss,2e-15)
    close('budget_width',row['width'],high[i]-low[i],2e-12)
    n=int(row['budget']);allocation=float(row['allocation'])
    assert int(row['n_plus'])==int(n*allocation)
    assert int(row['n_plus'])+int(row['n_minus'])==n
    checks['budget_rows']+=1

# Fresh replay of all 216,000 binomial counts from declared streams.
offset=0
for seed in (2026101001,2026101002,2026101003):
    random=np.random.default_rng(seed)
    for world,p,pm in [('A',29/60,19/30),('B',19/60,7/15),('C',.35,.5)]:
        for size in (200,800,3200,12800):
            for numerator,denominator in ((1,3),(1,2),(2,3)):
                n1=int(size*(numerator/denominator));n2=size-n1
                expected1=random.binomial(n1,p,size=1000);expected2=random.binomial(n2,pm,size=1000)
                subset=budget[offset:offset+1000]
                assert all(r['seed']==str(seed) and r['world']==world and int(r['budget'])==size for r in subset)
                assert np.array_equal(expected1,[int(r['k_plus']) for r in subset])
                assert np.array_equal(expected2,[int(r['k_minus']) for r in subset])
                offset+=1000
assert offset==len(budget)
checks['regenerated_binomial_counts']=2*offset

for name,data,groupkeys,seedkeys in [('joint',out,['world','recording','method'],['seed','world','recording','method']),('budget',budget,['world','budget','allocation'],['seed','world','budget','allocation'])]:
    summary=json.loads((NEW/(name+'_summary.json')).read_text(encoding='utf-8'))
    summary_check(data,summary['summary'],groupkeys)
    summary_check(data,summary['per_seed'],seedkeys)

# Independent randomized forward-law identities; full mechanism generation then inverse.
random=np.random.default_rng(7011010)
for _ in range(2000):
    b=random.uniform(0,.3);s=random.uniform(b+.5,1);p=random.uniform(0,1);lam=random.uniform(1,8)
    q=p*s+(1-p)*b
    joint_positive=p*s
    weighted_diseased=p+(lam-1)*joint_positive
    selected=weighted_diseased/(1+(lam-1)*q)
    restored=inverse_odds(selected,lam,s,b)
    close('forward_inverse_identity',restored,p,3e-15)
    checks['forward_inverse_identity']+=1

grid_delta=max(max(float(r['128_lower'])-float(r['64_lower']),float(r['64_upper'])-float(r['128_upper'])) for r in gridrows)
report=dict(status='PASS_DETERMINISTIC',review_independence='same-family',acceptance_status='provisional',
            implementation_independence='Separate code; no imports from analysis scripts; shared NumPy/SciPy numerical libraries',
            checks=dict(checks),maximum_absolute_discrepancies=dict(max_error),
            maximum_64_to_128_endpoint_change=grid_delta,elapsed_seconds=time.time()-started,
            scope='E3/E4 numerical and assumption audit; not clinical validity, general optimality, external peer review, or final manuscript readiness')
(ROOT/'verification/INDEPENDENT_EXTENSION_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)

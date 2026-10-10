"""Exact population calculation; protocol: GROUP_CHANNEL_PROTOCOL.md."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import csv, json
import numpy as np
from scipy.optimize import linprog

OUT=Path(__file__).resolve().parent
Q=(F(31,100),F(2,5)); P0=(F(2,5),F(11,20))
DP=list(map(F,['0','.01','.02','.03','.05','.075','.10','.15']))
DC=list(map(F,['0','.005','.01','.02','.03','.04','.05']))
WORLDS={'A':(F(29,60),F(19,30)), 'B':(F(19,60),F(7,15)), 'C':(F(7,20),F(1,2))}

def polygon(dp,dc,groups,audits):
    rows=[(-1,0,-F(1,2)),(0,-1,-(F(7,100)-dc)),(0,1,F(7,100)+dc),
          (-1,-1,-(F(67,100)-dc)),(1,1,F(67,100)+dc),(0,-1,0),(0,1,1),(-1,-1,0),(1,1,1)]
    for g in groups:
        lo=max(F(0),P0[g]-dp); hi=min(F(1),P0[g]+dp)
        if g in audits: lo=max(lo,audits[g]); hi=min(hi,audits[g])
        rows.extend([(lo,1,Q[g]),(-hi,-1,-Q[g])])
    rows=[tuple(map(F,row)) for row in rows]
    vertices=set()
    for (a,b,c),(d,e,f) in combinations(rows,2):
        det=a*e-b*d
        if det:
            x=(c*e-b*f)/det; y=(a*f-c*d)/det
            if all(h*x+j*y<=k for h,j,k in rows): vertices.add((x,y))
    return sorted(vertices)

def witness(channels):
    pp=[(Q[g]-b)/a for g,(a,b) in enumerate(channels)]
    return {'risks':[str(x) for x in pp], 'sensitivity':[str(a+b) for a,b in channels],
            'specificity':[str(1-b) for a,b in channels]}

def exact(dp,dc,shared,audits):
    if shared:
        vv=polygon(dp,dc,(0,1),audits)
        candidates=[((10*(Q[0]-b)-5*(Q[1]-b))/a-1)/12 for a,b in vv]
        ww=[witness((v,v)) for v in vv]
    else:
        vv=[polygon(dp,dc,(g,),audits) for g in (0,1)]
        choices=[]
        for g in (0,1):
            if not vv[g]: return None
            ordered=sorted(vv[g],key=lambda v:(Q[g]-v[1])/v[0])
            choices.append((ordered[0],ordered[-1]))
        ww=[witness((choices[0][0],choices[1][1])),witness((choices[0][1],choices[1][0]))]
        candidates=[(10*F(w['risks'][0])-5*F(w['risks'][1])-1)/12 for w in ww]
    if not candidates: return None
    low=min(range(len(candidates)),key=candidates.__getitem__)
    high=max(range(len(candidates)),key=candidates.__getitem__)
    return candidates[low],candidates[high],ww[low],ww[high]

def lp(dp,dc,shared,audits):
    # An independent transformed formulation; variables are u,v,p+,(u-,v-),p-.
    n=4 if shared else 6; rows=[]; rhs=[]
    def add(items,b):
        r=np.zeros(n)
        for i,x in items.items(): r[i]=float(x)
        rows.append(r); rhs.append(float(b))
    pi=(2,3) if shared else (2,5)
    for g in range(1 if shared else 2):
        ui,vi=(0,1) if shared or g==0 else (3,4)
        sl=F(67,100)-dc; su=F(67,100)+dc
        bl=F(7,100)-dc; bu=F(7,100)+dc
        add({ui:sl,vi:-1},1); add({ui:-su,vi:1},-1)
        add({ui:bl,vi:-1},0); add({ui:-bu,vi:1},0)
    equal=[]; erhs=[]
    bounds=[(None,None)]*n
    for g in (0,1):
        ui,vi=(0,1) if shared or g==0 else (3,4)
        bounds[ui]=(1,2); bounds[vi]=(0,None)
        lo=max(F(0),P0[g]-dp); hi=min(F(1),P0[g]+dp)
        if g in audits: lo=max(lo,audits[g]); hi=min(hi,audits[g])
        if lo>hi: return None
        bounds[pi[g]]=(float(lo),float(hi))
        r=np.zeros(n); r[ui]=-float(Q[g]); r[vi]=1; r[pi[g]]=1
        equal.append(r); erhs.append(0)
    c=np.zeros(n); c[pi[0]]=10/12; c[pi[1]]=-5/12
    r=[linprog(sign*c,A_ub=rows,b_ub=rhs,A_eq=equal,b_eq=erhs,bounds=bounds,method='highs') for sign in (1,-1)]
    if r[0].status==2: return None
    assert all(x.success for x in r),[x.message for x in r]
    return r[0].fun-1/12,-r[1].fun-1/12

def run():
    records=[]; max_error=0
    def cell(dp,dc,shared,audits,label):
        nonlocal max_error
        result=exact(dp,dc,shared,audits); numeric=lp(dp,dc,shared,audits)
        assert (result is None)==(numeric is None)
        row={'risk_cap':str(dp),'channel_cap':str(dc),'channel':'shared' if shared else 'group_specific','audit':label}
        if result is None: row['status']='infeasible'
        else:
            lo,hi,lw,hw=result
            error=max(abs(float(lo)-numeric[0]),abs(float(hi)-numeric[1])); max_error=max(max_error,error)
            assert error<1e-9,(row,result,numeric)
            for w in (lw,hw):
                p=list(map(F,w['risks'])); se=list(map(F,w['sensitivity'])); sp=list(map(F,w['specificity']))
                for g in (0,1):
                    assert se[g]*p[g]+(1-sp[g])*(1-p[g])==Q[g]
                    assert abs(p[g]-P0[g])<=dp and abs(se[g]-F(67,100))<=dc and abs(sp[g]-F(93,100))<=dc
                    assert se[g]+sp[g]-1>=F(1,2)
                    if g in audits: assert p[g]==audits[g]
                if shared: assert se[0]==se[1] and sp[0]==sp[1]
            row.update(lower=str(lo),upper=str(hi),lower_decimal=float(lo),upper_decimal=float(hi),width=str(hi-lo),
                       status='positive' if lo>0 else 'negative' if hi<0 else 'zero_compatible',lower_witness=lw,upper_witness=hw)
        records.append(row)
        return result
    grid={}
    for dp in DP:
        for dc in DC:
            a=cell(dp,dc,True,{},'none'); b=cell(dp,dc,False,{},'none')
            assert b[0]<=a[0]<=a[1]<=b[1]
            grid[dp,dc,True]=a; grid[dp,dc,False]=b
    for (dp,dc,s),v in grid.items():
        for (dp2,dc2,s2),w in grid.items():
            if s==s2 and dp2>=dp and dc2>=dc: assert w[0]<=v[0]<=v[1]<=w[1]
    for name,risks in WORLDS.items():
        for selected in ((0,),(1,),(0,1)):
            aa={g:risks[g] for g in selected}
            vals=[]
            for shared in (True,False):
                v=cell(F(1,10),F(1,20),shared,aa,name+'_'+''.join('plus' if g==0 else 'minus' for g in selected))
                base=grid[F(1,10),F(1,20),shared]
                assert base[0]<=v[0]<=v[1]<=base[1]
                vals.append(v)
                if len(selected)==2: assert v[0]==v[1]==(10*risks[0]-5*risks[1]-1)/12
            assert vals[1][0]<=vals[0][0]<=vals[0][1]<=vals[1][1]
    result={'protocol':'GROUP_CHANNEL_PROTOCOL.md','kind':'exact_synthetic_population','rows':records,
            'checks':{'cells':len(records),'lp_max_absolute_error':max_error,'exact_witness_and_nesting':'PASS'}}
    (OUT/'group_channels_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    keys=['risk_cap','channel_cap','channel','audit','lower','upper','lower_decimal','upper_decimal','width','status']
    with (OUT/'group_channels_results.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore'); writer.writeheader(); writer.writerows(records)
    print(json.dumps(result['checks']))
    for r in records:
        if r['risk_cap']=='1/10' and r['channel_cap']=='1/20': print({k:r[k] for k in keys if k in r})

if __name__=='__main__': run()

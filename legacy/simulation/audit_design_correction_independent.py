"""Independent all-row numerical audit of the post hoc design correction.

No production simulation, projection, or correction functions are imported.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import csv
import hashlib
import json
import math
from pathlib import Path

from scipy.optimize import linprog
from scipy.special import betaincinv

D = Path(__file__).resolve().parent
N = 10_000
A = .05 / 7
TOL = 1e-10
P = ("source_gold_plus", "source_gold_minus", "source_se", "source_sp",
     "target_proxy_plus", "target_proxy_minus", "target_audit_recorded")
W = {
    "B_negative": (F(19,60), F(7,15), F(18,25), F(22,25), F(-1,72),
                   (F(-1,72),F(-1,72))),
    "C_interior_tie": (F(7,20), F(1,2), F(7,10), F(9,10), F(0),
                       (F(-1,152),F(1,328))),
}
KEYS = {(w,l,u,v) for w in W for l in (1.,4.)
        for u,v in ((1.,1.),(.8,.95))}

def sha(name):
    return hashlib.sha256((D/name).read_bytes()).hexdigest().upper()

def close(a,b,tol=TOL):
    return math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tol

@lru_cache(maxsize=None)
def cp(x,n,alpha=A):
    assert 0 <= x <= n
    if n == 0: return (0.,1.)
    return (
        0. if x == 0 else float(betaincinv(x,n-x+1,alpha/2)),
        1. if x == n else float(betaincinv(x+1,n-x,1-alpha/2)),
    )

def key(row):
    return row["world"],float(row["lambda_proxy_positive"]),float(row["ref_se"]),float(row["ref_sp"])

def truth(k):
    p,m,se,sp,d,_ = W[k[0]]
    lam,u,v = map(lambda x:F(str(x)),k[1:])
    q = F(31,100)
    r = p*(1+(lam-1)*se)/(1+(lam-1)*q)
    a = (u+v-1)*r+1-v
    return r,a

def inversion(a,q,se,lam,u,v):
    al,ah=max(a[0],1-v),min(a[1],u)
    if al > ah: return None
    sl,su=max(0.,se[0]-.05),min(1.,se[1]+.05)
    gain=u+v-1
    k=lam-1
    lo=max(0.,((al+v-1)/gain)*(1+k*q[0])/(1+k*su))
    hi=min(1.,((ah+v-1)/gain)*(1+k*q[1])/(1+k*sl))
    if lo > hi + 1e-12: return None
    return lo,hi

def lp(box,audit):
    """p+,p-,h,z with h=1/(Se+Sp-1), z=(1-Sp)h; qh=p+z."""
    gp,gm,se,sp,qp,qm,_=box
    sl,su=max(0.,se[0]-.05),min(1.,se[1]+.05)
    cl,cu=max(0.,sp[0]-.05),min(1.,sp[1]+.05)
    pl,ph=max(0.,gp[0]-.1,audit[0]),min(1.,gp[1]+.1,audit[1])
    ml,mh=max(0.,gm[0]-.1),min(1.,gm[1]+.1)
    if pl>ph or ml>mh: return None
    aa=[
      [0.,0.,sl,-1.],[0.,0.,-su,1.],
      [0.,0.,1-cu,-1.],[0.,0.,cl-1,1.],
      [-1.,0.,qp[0],-1.],[1.,0.,-qp[1],1.],
      [0.,-1.,qm[0],-1.],[0.,1.,-qm[1],1.],
    ]
    bb=[1.,-1.,0.,0.,0.,0.,0.,0.]
    c=[5/6,-5/12,0.,0.]
    bounds=[(pl,ph),(ml,mh),(1.,2.),(None,None)]
    low=linprog(c,A_ub=aa,b_ub=bb,bounds=bounds,method="highs")
    if low.status==2: return None
    assert low.success,(low.status,low.message)
    high=linprog([-z for z in c],A_ub=aa,b_ub=bb,bounds=bounds,method="highs")
    assert high.success,(high.status,high.message)
    return float(low.fun-1/12),float(-high.fun-1/12)

def compare_lp(row,projected,corrected):
    if projected is None:
        assert int(row["infeasible"])==1
        for f in ("lower","upper","width"):
            if corrected: assert row[f]=="",(f,row[f])
            else: assert math.isnan(float(row[f])),(f,row[f])
        return 0.
    assert int(row["infeasible"])==0
    d=max(abs(float(row["lower"])-projected[0]),
          abs(float(row["upper"])-projected[1]))
    assert d<=TOL,d
    assert close(float(row["width"]),projected[1]-projected[0])
    return d

def main():
    names=("DESIGN_CORRECTION_PROTOCOL_V1.md",
           "PRE_RESULT_DESIGN_CORRECTION_REVIEW.md",
           "run_design_correction_v1.py","run_supplement_v4.py",
           "supplement_v3_raw_main_10000.csv",
           "supplement_v4_raw_main_10000.csv",
           "design_correction_v1_raw.csv",
           "design_correction_v1_summary.json")
    hashes={name:sha(name) for name in names}
    summary=json.loads((D/"design_correction_v1_summary.json").read_text(encoding="utf-8"))
    assert summary["n_total"]==80_000 and summary["n_per_cell"]==N
    assert summary["status"]=="POST_HOC_EXPLORATORY_MECHANISM_KNOWN"
    assert summary["protocol_sha256"]==hashes["DESIGN_CORRECTION_PROTOCOL_V1.md"]
    assert summary["script_sha256"]==hashes["run_design_correction_v1.py"]
    assert summary["input_sha256"]==hashes["supplement_v4_raw_main_10000.csv"]
    assert summary["output_sha256"]==hashes["design_correction_v1_raw.csv"]
    cells={(c["world"],c["lambda_proxy_positive"],c["ref_se"],c["ref_sp"]):c
           for c in summary["cells"]}
    assert set(cells)==KEYS
    seen=defaultdict(set)
    counts=defaultdict(lambda:defaultdict(int))
    widths=defaultdict(float)
    widths_n=defaultdict(int)
    maxd=defaultdict(float)
    total=0
    with (D/"supplement_v4_raw_main_10000.csv").open(newline="",encoding="utf-8") as fa, \
         (D/"design_correction_v1_raw.csv").open(newline="",encoding="utf-8") as fb, \
         (D/"supplement_v3_raw_main_10000.csv").open(newline="",encoding="utf-8") as fc:
        for old,new,prior in zip(csv.DictReader(fa),csv.DictReader(fb),
                                 csv.DictReader(fc),strict=True):
            total+=1
            k=key(old)
            assert k in KEYS and key(new)==k
            assert all(old[f]==prior[f] for f in prior),("V3/V4 mismatch",total)
            rep=int(old["rep"])
            assert rep==int(new["rep"]) and 0<=rep<N and rep not in seen[k]
            seen[k].add(rep)
            p,m,se,sp,d,control_set=W[k[0]]
            assert se*p+(1-sp)*(1-p)==F(31,100)
            assert se*m+(1-sp)*(1-m)==F(2,5)
            selected,recorded=truth(k)
            for f,v in (("audit_selected_true_risk",selected),
                        ("audit_recorded_positive_rate",recorded)):
                dv=abs(float(old[f])-float(v))
                maxd["generator_rates"]=max(maxd["generator_rates"],dv)
                assert dv<=TOL,(total,f,dv)
            for row in (old,new):
                dv=abs(float(row["true_delta"])-float(d))
                maxd["true_delta"]=max(maxd["true_delta"],dv)
                assert dv<=TOL
            nplus=int(old["source_gold_plus_success"])
            nminus=int(old["source_gold_minus_success"])
            denominators=(1000,1000,nplus+nminus,2000-nplus-nminus,1000,1000,800)
            box=[]
            for prim,n_expect in zip(P,denominators):
                x=int(old[prim+"_success"])
                n=int(old[prim+"_denominator"])
                assert n==n_expect and 0<=x<=n,(total,prim,x,n)
                calc=cp(x,n)
                saved=(float(old[prim+"_lower"]),float(old[prim+"_upper"]))
                dv=max(abs(calc[0]-saved[0]),abs(calc[1]-saved[1]))
                maxd["seven_cp"]=max(maxd["seven_cp"],dv)
                assert dv<=1e-12,(total,prim,dv)
                box.append(calc)
            naive=lp(box,box[6])
            maxd["naive_lp"]=max(maxd["naive_lp"],compare_lp(old,naive,False))
            correction=inversion(box[6],box[4],box[2],*k[1:])
            if correction is None:
                assert new["corrected_p_lower"]=="" and new["corrected_p_upper"]==""
                revised=None
            else:
                dv=max(abs(float(new["corrected_p_lower"])-correction[0]),
                       abs(float(new["corrected_p_upper"])-correction[1]))
                maxd["inverse"]=max(maxd["inverse"],dv)
                assert dv<=TOL and correction[0]<=correction[1]+TOL
                revised=lp(box,correction)
            maxd["corrected_lp"]=max(maxd["corrected_lp"],compare_lp(new,revised,True))
            control=k[1:]==(1.,1.,1.)
            if control:
                assert (naive is None)==(revised is None),("identity feasibility",total)
                if naive is not None:
                    dv=max(abs(naive[0]-revised[0]),abs(naive[1]-revised[1]))
                    maxd["identity_lp"]=max(maxd["identity_lp"],dv)
                    assert dv<=TOL
            if revised is None:
                flags={"true_covered":0,"correct_sign":0,"wrong_sign":0,
                       "indeterminate":0,"infeasible":1}
                set_covered=0
            else:
                lo,hi=revised
                pos=lo>TOL
                neg=hi<-TOL
                assert not (pos and neg)
                correct=(pos and d>0) or (neg and d<0)
                wrong=(pos and d<=0) or (neg and d>=0)
                flags={"true_covered":int(lo-TOL<=float(d)<=hi+TOL),
                       "correct_sign":int(correct),"wrong_sign":int(wrong),
                       "indeterminate":int(not(correct or wrong)),"infeasible":0}
                set_covered=int(lo-TOL<=float(control_set[0]) and
                                float(control_set[1])<=hi+TOL) if control else 0
                widths[k]+=hi-lo
                widths_n[k]+=1
            assert sum(flags[x] for x in ("correct_sign","wrong_sign",
                                           "indeterminate","infeasible"))==1
            for field,val in flags.items():
                assert int(new[field])==val,(total,field,val,new[field])
                counts[k][field]+=val
            if control:
                assert int(new["control_set_covered"])==set_covered
                counts[k]["control_set_covered"]+=set_covered
            else:
                assert new["control_set_covered"]==""
            naive_covered=int(naive is not None and
                              naive[0]-1e-9<=float(d)<=naive[1]+1e-9)
            naive_wrong=int(naive is not None and
                            ((naive[0]>1e-12 and d<=0) or
                             (naive[1]<-1e-12 and d>=0)))
            assert int(old["true_covered"])==naive_covered
            assert int(old["wrong_sign"])==naive_wrong
            counts[k]["naive_true_covered"]+=naive_covered
            counts[k]["naive_wrong_sign"]+=naive_wrong
            if total%10_000==0: print("AUDITED",total,flush=True)
    assert total==8*N and set(seen)==KEYS and all(v==set(range(N)) for v in seen.values())
    for k in sorted(KEYS):
        metrics=cells[k]["metrics"]
        for field in ("true_covered","correct_sign","wrong_sign",
                      "indeterminate","infeasible","naive_true_covered",
                      "naive_wrong_sign")+ (("control_set_covered",) if k[1:]==(1.,1.,1.) else ()):
            x=counts[k][field]
            assert metrics[field]["count"]==x,(k,field)
            assert close(metrics[field]["proportion"],x/N,1e-14)
            ci=cp(x,N,.05)
            assert all(close(y,z,1e-12) for y,z in
                       zip(metrics[field]["mc_95_interval"],ci)),(k,field,"CI")
        expected_width=widths[k]/widths_n[k] if widths_n[k] else None
        if expected_width is None: assert metrics["mean_width_feasible_only"] is None
        else: assert close(metrics["mean_width_feasible_only"],expected_width,1e-12)
    result={
      "status":"PASS_DETERMINISTIC",
      "review_independence":"fresh_same_family_provisional",
      "n_total":total,"n_cells":8,"n_per_cell":N,
      "verified":["generator truth and all eight mechanism probabilities",
                  "V3/V4 source row identity","seven CP intervals per row",
                  "naive and corrected LPs per row","correction inversion per row",
                  "all flags and summary counts/Monte Carlo intervals",
                  "control identity including feasibility","unique complete replicate IDs"],
      "max_abs_discrepancy":dict(maxd),
      "sha256":hashes,
      "generator_truth":[],
      "cells":[],
      "limitations":["synthetic mechanism-known scenario only",
                     "fresh same-family audit is provisional",
                     "mechanism-specific sharp population sets outside controls were not derived",
                     "no real-world clinical performance or mechanism-parameter uncertainty evaluated"],
      "blocking_interpretation":["do not claim whole-set coverage outside controls",
                                 "do not call the 95% Monte Carlo intervals patient-level uncertainty",
                                 "post hoc exploratory analysis; author/editorial signoff separate"],
    }
    for k in sorted(KEYS):
        r,a=truth(k)
        p,m,se,sp,d,_=W[k[0]]
        result["generator_truth"].append({
            "world":k[0],"lambda":k[1],"ref_se":k[2],"ref_sp":k[3],
            "p_plus":str(p),"p_minus":str(m),"se":str(se),"sp":str(sp),
            "q_plus":"31/100","q_minus":"2/5","selected_risk":str(r),
            "recorded_probability":str(a),"true_delta":str(d)})
        result["cells"].append({
            "world":k[0],"lambda":k[1],"ref_se":k[2],"ref_sp":k[3],
            "n":len(seen[k]),"counts":dict(counts[k]),
            "mean_width_feasible_only":widths[k]/widths_n[k] if widths_n[k] else None})
    (D/"DESIGN_CORRECTION_INDEPENDENT_AUDIT.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("PASS_DETERMINISTIC",json.dumps(dict(maxd)),flush=True)

if __name__=="__main__":
    main()

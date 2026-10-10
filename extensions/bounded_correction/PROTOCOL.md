# Bounded audit-mechanism correction: frozen extension protocol

Frozen 30 September 2026 before extension code or outcomes. Same-family provisional design review. This new synthetic extension is motivated by already known original stress-test results; it is not a prospective clinical validation or a retroactive addition to the original ADEMP lock. Original companion files remain unchanged.

## Claims and falsifiable aims

Primary claim: the interpretation of an uncertainty interval depends on whether selection and reference-recording assumptions are credible. A correction based on a true-containing mechanism class can remove unsupported precision, with possible loss of sign certification. Supporting claim: mechanism uncertainty cannot be inferred away from the same seven aggregate counts.

Aim 1: verify that the proposed conservative inverse rectangle contains the generating target risk whenever the seven simultaneous primitive intervals cover and the mechanism box contains the generating parameters. Any counterexample, missing row, or failed LP certificate is a failure. Aim 2: quantify true-point coverage, wrong nonzero-sign certification, infeasibility and width under known versus bounded versus false mechanism assumptions. No minimum superiority, narrower-width or power result is required. Negative or null findings must be retained.

## Frozen design

Use original B (p+ = 19/60, p- = 7/15, target Se = .72, Sp = .88, Delta = -1/72) and C (p+ = .35, p- = .50, target Se = .70, Sp = .90, Delta = 0). Both have target proxy probabilities (.31, .40). Fixed weights (2/3, 1/3), threshold .20; source risks (.40, .55), source proxy Se .67, Sp .93; risk drift .10, channel drift .05, target Youden at least .50 and a common target proxy channel remain unchanged.

Four generating cells: B/C crossed with accurate audit recording (u,v)=(1,1) or fallible recording (.80,.95); selection multiplier lambda=4 in all cells. Sample 1,000 source paired reference/proxy outcomes per group, 1,000 target proxy outcomes per group and a separate 800-observation selected G+ audit. Draw the seven binomial sufficient statistics directly; do not generate the discarded historical V3 audit draw. Store all counts and denominators. Selection is independent weighted sampling from the population, not sampling without replacement from a finite cohort.

Three independent NumPy PCG64 roots: 2026093001, 2026093002, 2026093003. For each root spawn four child streams in B-accurate, B-fallible, C-accurate, C-fallible order. Each cell/seed has 1,000 repetitions: 12,000 new evidence records, 3,000 per generating cell. Compare four methods on identical records (paired comparisons):

| Analysis | Lambda assumption | Audit sensitivity u | Audit specificity v |
|---|---|---|---|
| Naive | 1 | 1 | 1 |
| Oracle-known | 4 | generating value | generating value |
| True-containing bounded | [3,5] | [.95,1] if accurate; [.75,.85] if fallible | [.95,1] if accurate; [.90,1] if fallible |
| Wrong-selection-bound control | [1,1.25] | generating value | generating value |

These numerical ranges are diagnostic stipulations, not clinical estimates or evidence that the analyst could identify the channel. The wrong-bound control isolates an erroneous selection restriction; its interval need not be valid, and no particular failure rate is required.

## Mechanism assumptions and mapping

Let Y be reference truth, Z the routine proxy, S audit inclusion and R the audit record. P(S=1|Y,Z)/P(S=1|Y,Z=0) is lambda when Z=1: S is independent of Y conditional on Z within G+. The constant reference channel obeys P(R=1|Y=1,S=1,Z)=u and P(R=0|Y=0,S=1,Z)=v. Thus recording is independent of Z conditional on Y among selected subjects, with the same u,v across selected proxy strata. Differential review error or further outcome-dependent selection is outside this extension's guarantee.

Write q=P(Z=1), s=P(Z=1|Y=1), p=P(Y=1), a=P(R=1|S=1). Then r=P(Y=1|S=1)=p[1+(lambda-1)s]/[1+(lambda-1)q] and a=(u+v-1)r+1-v. Mechanism parameters lambda,u,v are not identified by the available aggregate q,a and source counts. In practice their ranges need external design records, validation or explicit sensitivity elicitation. If estimated with confidence intervals, allocate additional simultaneous error probability; the present boxes are deterministic restrictions and consume no extra statistical alpha.

For each of seven primitives use two-sided Clopper-Pearson at error .05/7. Random source sensitivity/specificity denominators are handled conditionally; a zero denominator gives [0,1]. Independence among primitive intervals is not required for the union bound. Expand the source Se interval by .05 and clip to [0,1] to bound target s.

For box u in [ul,uh], v in [vl,vh] with ul+vl>1, and audit interval [al,ah], reject if al>uh or ah<1-vh. A conservative selected-risk interval is

rL=max(0,(al+vl-1)/(uh+vl-1)), rU=min(1,(ah+vh-1)/(ul+vh-1)).

The inversion increases with a and v and decreases with u on the physically feasible domain. For lambda in [ll,lh], evaluate H=(1+(lambda-1)q)/(1+(lambda-1)s) at both lambda endpoints. Hmin uses qL,sU; Hmax uses qU,sL. This endpoint rule is valid because its derivative in lambda has constant sign q-s for fixed q,s. Then pL=max(0,rL Hmin), pU=min(1,rU Hmax). Reject an empty interval. Intersect it with the unchanged common-channel projection constraints. This is an outer envelope, not a sharp bounded-mechanism identified set; it can discard joint compatibility and be conservative.

The oracle is a singleton mechanism box. The naive singleton reduces exactly to the untransformed audit interval. The bounded true-containing box must contain the oracle inverse interval and, if oracle projection is feasible, the bounded projection must be feasible and contain its endpoints within 1e-9. This rowwise nesting is required independently of favorable outcomes.

## Reporting and validation

Primary metrics: true-point coverage and wrong nonzero-sign certification. Secondary: correct-sign certification for B, indeterminate and infeasible fractions, feasible-only width, paired bounded-minus-oracle width where both feasible. For C, any nonzero sign is wrong; abstention is appropriate and is not labeled correct-sign certification. All 3,000 rows per generating cell remain in frequency denominators; infeasible rows count as noncoverage and never force a sign.

Report each seed separately and pooled cell counts, proportions, MCSE and exact two-sided 95% binomial Monte Carlo intervals; for zero wrong calls also give a one-sided 95% upper limit. These intervals quantify Monte Carlo error, not patient-level uncertainty and not simultaneous multiplicity-adjusted intervals across reported metrics. At n=3,000, worst-case MCSE is 0.00913, MCSE at .95 is 0.00398, and a zero-event one-sided 95% upper limit is about .000998. No optional stopping, sample-size increase, scenario alteration or favorable-seed selection is allowed after results.

Require deterministic boundary and nesting tests before simulations; record protocol/script/input/output hashes and package versions. Save every primitive count, interval, mechanism metadata and method endpoint/flag. An independent audit must reconstruct all primitive intervals, corrections and projections without importing the implementation, verify unique complete cell/seed/rep/method keys, and validate all summaries. Direct deterministic population mappings are labeled non-sharp outer envelopes. Do not report mechanism-specific whole-set simulation coverage because those sharp sets are not derived here.

Nominal coverage is a theorem conditional on correct sampling and true-containing restrictions; simulation is a selected-mechanism implementation check. Reject nominal empirical validity only if a multiplicity-adjusted upper confidence bound for a valid method's coverage is below .95, or if a deterministic containment invariant fails. Above-.95 finite-cell performance is not proof of uniform validity. Methods with false restrictions carry no nominal guarantee.

## Scope and execution

CPU only; approximately 48,000 method envelopes and 96,000 small LPs plus deterministic checks. A 10-row-per-cell/seed engineering pilot may test runtime/schema; keep it separate and do not use it to revise scientific choices. Run all frozen cells, then independent audit. No patient data, ML fitting, efficiency ranking, optimized audit allocation or clinical calibration is introduced.

# Locked supplemental analysis protocol

Date: 10 October 2026, before executing this extension. This is a revision-stage, hypothesis-informed extension after seeing earlier results, not a prospective registration of the original study.

## E3 Joint audit information

Aim: distinguish loss of information in separate interval inversion from genuine identification ambiguity. Reuse ALL 12000 evidence records from the 30 September bounded-correction experiment, without resampling or dropping records. Compare old oracle and bounded outer envelopes with a joint, interval-cell outer projection using the same seven Bonferroni-adjusted Clopper-Pearson primitives and the same mechanism boxes. B/C worlds crossed with accurate/fallible reference, three seeds of 1000 per cell. No new clinical observations.

For target sensitivity s, false-positive probability b, J=s-b>=0.5, q=b+Jp. If selection odds ratio is lambda>=1 and selected true reference risk is r, then p=r[1+(lambda-1)b]/{1+(lambda-1)[(1-r)s+rb]}. Invert recording error to a feasible r interval; retain this equation jointly with proxy and source bounds through conservative cell bounds over (s,b). Use a 64 by 64 partition for every record, intersecting with the old envelope because both contain the same feasible set. Report as a conservative outer envelope, never a sharp finite-sample set. Use 1e-10 outward padding. Empty projections remain failures/infeasible, not omitted from coverage denominators.

Performance: true-point coverage, correct/wrong sign certification, indeterminate and infeasible fractions; mean widths (denominator explicit), paired width differences; Monte Carlo uncertainty on proportions and mean differences. Reusing samples gives paired analysis, not additional independent replications. Verify true-point containment for all records whose primitives jointly cover and bounded-versus-oracle nesting.

Numerical sensitivity: first 20 records in each seed/cell (240 records), selected by file order independent of outcomes, compare 32/64/128 partitions; retain all differences. Population cells use 64/128/256, and exact oracle B/accurate singleton -1/72 is an analytic anchor. Mesh intervals are outer guarantees up to floating-point implementation tolerance; refinement agreement does not establish sharpness.

## E4 Equal incremental reference budget

Aim: compare declared audit allocations and distinguish uniform sign certification from a decision under a stated loss. Direct representative accurate audits in BOTH disagreement strata avoid source-to-target and shared-channel assumptions. Known weights (2/3,1/3), t=0.2, worlds A/B/C from the original model. The historical source/proxy information is not used by any compared procedure. This comparison concerns new reference ascertainments only; equal per-case costs and complete reference ascertainment are stipulated, and screening/recruitment costs are excluded.

Budgets N=200,800,3200,12800 total reference outcomes. Fixed plus allocations 1/3,1/2,2/3 (integer floor, rest to minus); no selection after viewing results, no claim of optimal allocation. 3000 repetitions per world/budget/allocation, split across seeds 2026101001, 2026101002, 2026101003 with 1000 per cell each. Independent binomial sampling per stratum. Two exact binomial intervals each at 97.5% produce a >=95% simultaneous rectangle and the linear net-benefit contrast envelope.

Report coverage, correct/wrong sign, mean width, minimax-regret choice based on the entire envelope (choose rule 1 if max(0,-L)<=max(0,U), otherwise rule 2), empirical true net-benefit loss relative to the better fixed rule, and Monte Carlo errors. At an exact tie any strict sign is false; either decision has zero regret. Do not equate a nonsignificant result with a prohibition on choosing. These empirical losses are conditional operating characteristics, not expected clinical value, acquisition net value or prospective policy validation.

## Checks and output

Archive raw counts, endpoints, decisions and all cells, protocol, code, environment and numerical checks. Independent reviewer checks derivation and separate-code reconstruction. Keep all adverse findings. Original source/audit assumptions and the reference-outcome limitations are unchanged. New scientific claims must be limited to these designs.

## Terminology correction after independent review
The E3 phrase selection odds ratio above means the ratio of inclusion probabilities for proxy-positive versus proxy-negative cases (selection-weight multiplier). It is not the odds ratio of inclusion. This wording correction changes no equation, simulation, mechanism value or result. Original locked text and its execution hash are retained.


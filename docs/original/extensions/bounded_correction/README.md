# Reproducible bounded audit-mechanism extension

This directory contains a synthetic extension to the original manuscript's audit example. It does not contain patient data. The protocol was fixed before the extension implementation and results; it is supplementary new work motivated by previously observed original results. No result is represented as part of the original prespecified experiment.

## Files

- `PROTOCOL.md` and its dated copy: immutable scientific design; SHA-256 `9E734121754BDAECCF2B67DE18A8E71F9FD64530BF59C39EE52B8A40B39AAFFD`.
- `run_bounded_correction.py`: generator, conservative correction and paired analysis; imports the archived original LP unchanged.
- `primitive_counts_main.csv`: 12,000 independently generated evidence records, 1,000 per seed and generating cell.
- `method_results_main.csv`: 48,000 method records; join to evidence on `seed,cell,rep` and distinguish analysis with `method`.
- `summary_main.json`: frozen analysis output with pooled and per-seed results, exact Monte Carlo intervals, pairwise width differences, source/data hashes and environment versions. Its status intentionally remains `COMPLETE_PENDING_INDEPENDENT_AUDIT`; subsequent validation is recorded separately rather than rewriting frozen results.
- `audit_bounded_correction.py`: separate reconstruction without importing simulation/correction code, including physical-corner inversion and reordered LP constraints.
- `AUDIT.json`: complete row-level audit result and hashes, once the audit has finished.
- `*_pilot.csv`, `summary_pilot.json`: engineering pilot only. These are excluded from the main counts and publication summaries. They use initial stream subsequences and are not additional independent repetitions.

## Reproduce

Use Python 3.12 with NumPy and SciPy. The actual run used Python 3.12.14, NumPy 2.3.5 and SciPy 1.18.1; this is a new run and does not claim random-stream identity with the original environment. From this directory:

```text
python run_bounded_correction.py --reps 10 --label pilot
python run_bounded_correction.py --reps 1000 --label main
python audit_bounded_correction.py
```

The analysis locates its unchanged LP dependency at companion-relative `legacy/simulation/run_study.py`. The packaged script has a path-only adaptation documented by packaging_manifest.json; its original execution copy is preserved separately. The audit depends only on the files within this directory and NumPy/SciPy. Run in a copy to preserve archived outputs; scripts intentionally write the named run files.

## Data dictionary and denominators

Evidence keys: `seed` is one of the three frozen PCG64 roots; `cell` 0/1 means B accurate/fallible recording and 2/3 means C accurate/fallible recording; `rep` runs 0–999. `world`, `recording` and the `*_true` columns record synthetic generating facts, not learned quantities. `selected_true_risk` is the audit population's reference risk; `audit_recorded_probability` includes audit reference error.

Each of seven primitive prefixes has `success`, `denominator`, `lower`, `upper`. Successes/denominators are the retained sufficient statistics. Bounds are two-sided Clopper–Pearson intervals at marginal error .05/7. Source proxy sensitivity/specificity denominators are conditional source reference-positive/negative totals. `all_primitives_cover` identifies the simulated event that all seven intervals contain their generating probabilities; it is available because the truth is stipulated and is not an observable diagnostic in a clinical study.

Method values are `naive`, `oracle`, `bounded` and `wrong_selection_bound`. `p_lower,p_upper` bound target G+ risk after inverse correction; `lower,upper,width` describe the projected net-benefit envelope. Empty endpoint/width fields indicate infeasibility, never zero width. All flags are 0/1. Exactly one of `correct_sign`, `wrong_sign`, `indeterminate`, `infeasible` is 1. For world C, every nonzero sign is wrong; abstention is `indeterminate`, not `correct_sign`.

All 3,000 records per generating cell remain in each method's pooled frequency denominator; the 1,000-per-seed summaries are also retained. Infeasible rows count as noncoverage. Mean width is explicitly feasible-only. Paired width differences are computed on identical evidence records where oracle and bounded projections are both feasible. Monte Carlo intervals describe simulation proportions, not individual risk and not the inferential coverage object.

## Interpretation boundary

The bounded boxes are stipulated sensitivity restrictions containing the chosen truth, not estimated clinical intervals. Selection depends only on the routine proxy within G+, and reference error is constant conditional on true outcome among selected subjects. The common target routine-proxy channel and original drift restrictions remain in force. The correction is a conservative outer projection, not a new sharp estimator. The audit supports computational reproducibility under these assumptions; it does not establish clinical efficacy, general validity under misspecification or cost efficiency.

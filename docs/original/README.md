# Computational companion: synthetic uncertainty and identification

Exact calculations and explicitly synthetic simulations. No patient information is included. This source-only release does not include the original completed-run data tables, publication packages or manuscript. Run in a copy because analysis commands overwrite their named outputs.

## Environments and entry points

Python 3.12.14 was used. Original and October simulations use NumPy2.5.3/SciPy1.18.1. September extensions used NumPy2.3.5/SciPy1.18.1; use their separate requirements_extensions.txt. Preserve these versions for seed-level stream reproduction. October plots use Matplotlib3.11.2, PyMuPDF1.28.2 and Arial; font substitutions can alter appearance. Install requirements_revision.txt for October numerical/figure reproduction.

Original exact analysis: `python legacy/analytic/check_exact.py` and `python legacy/analytic/check_threshold_corollary.py`. Legacy simulation commands and options are documented in legacy/README.txt and scripts beside the archived raw data. Their protocols explicitly distinguish original plans from subsequent extensions.

September exact channel extension: `python extensions/group_channels/group_channels.py`.
September complete correction simulation: `python extensions/bounded_correction/run_bounded_correction.py --reps 1000 --label main`.
September independent reconstruction: `python extensions/bounded_correction/audit_bounded_correction.py`.

October joint reanalysis plus equal-budget simulations: `python revision_20261010/experiments/run_extensions.py --part all`.
October independent all-row reconstruction: `python revision_20261010/verification/independent_extension_audit.py`.
Final publication figures: `python revision_20261010/figures/build_figures.py`. It regenerates1/6/7, reads archived September summary for6, and preserves verified vector panels2–5/S1/S2 with170mm PDF widths and600dpi TIFF/JPEG exports. Original panel plotting scripts and inputs remain in figures/ and legacy/.

## Data and chronology

E3 reuses all12000 September evidence records with two corrected methods (24000 method rows). E4 adds108000 simulations across3worlds,4budgets,3allocations and3seeds. They are not clinical samples. PROTOCOL_LOCKED_ORIGINAL.md is the executed pre-run protocol; PROTOCOL.md adds a terminology erratum: lambda is a ratio of inclusion probabilities, not a selection odds ratio. The original environment hash refers to the original protocol; a replay hashes the current corrected text. Neither the correction nor portable path changes alter the simulations.

October joint_results.csv identifies seed/world/recording/rep/method and reports bounds, coverage, sign/error indicators and paired width reduction. grid_sensitivity.csv reports32/64/128 intervals for a prespecified subset; joint_population.csv reports64/128/256 outer enclosures. budget_results.csv includes sample sizes, binomial counts, true contrast, allocation, bounds, coverage, sign indicators, minimax choice and true regret. Summary JSONs contain pooled and per-seed counts, exact binomial Monte Carlo intervals and mean standard errors. Legacy DATA_DICTIONARY.md documents earlier outcomes. Infeasible, tie, null and misspecified cases are retained.

Joint projection is a conservative interval-cell outer approximation, not a sharp solver or formal interval arithmetic. The budget comparison stipulates accurate reference ascertainment, equal per-case cost and known weights; it does not establish clinical effectiveness or an optimal study policy. Internal independent reconstructions share numerical libraries and automated model family; they are not external peer review.


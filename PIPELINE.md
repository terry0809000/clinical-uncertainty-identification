# Execution interfaces — Clinical prediction uncertainty and identification

This guide describes the preserved source interfaces by reading the code. No research module was imported, no command was run, and no performance or portability validation is implied. Use the original project instructions for the scientific run order.

## Data → analysis → outputs

**Input contract:** Analytical examples and explicitly synthetic simulations; no patient dataset required for the analytical checks.

1. Acquire inputs from their original providers and configure the paths defined in the source.
2. Install the recorded dependencies in an isolated environment.
3. Preserve grouping, seeds, development/evaluation separation and any protocol/hash gates.
4. Run acquisition/extraction and preprocessing before model development. Run evaluation only after the selection/freeze required by the original implementation.
5. Keep raw data, row-level derivatives, predictions and fitted models outside the public repository.

## Script reference

### `extensions/bounded_correction/audit_bounded_correction.py`

[Open source](extensions/bounded_correction/audit_bounded_correction.py) · Executable entry point

Separate all-row reconstruction: imports no simulation/correction implementation. Same-family algorithmic audit, not external independent authorship. Uses a separate reference-corner inversion and reordered LP, solved with highs-ds.

### `extensions/bounded_correction/run_bounded_correction.py`

[Open source](extensions/bounded_correction/run_bounded_correction.py) · Executable entry point

Frozen synthetic bounded audit correction; see PROTOCOL.md. CPU only.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--reps` | Yes | ; choices: (10,1000) |
| `--label` | Yes | ; choices: ('pilot','main') |

### `extensions/group_channels/group_channels.py`

[Open source](extensions/group_channels/group_channels.py) · Executable entry point

Exact population calculation; protocol: GROUP_CHANNEL_PROTOCOL.md.

### `figures/figure1_evidence_design.py`

[Open source](figures/figure1_evidence_design.py) · Supporting module

Editable graphical overview. Geometric icons; no generated patient imagery.

### `figures/figure2_identification.py`

[Open source](figures/figure2_identification.py) · Supporting module

Figure 2: exact analytic identified sets for observationally equivalent worlds. The worlds are synthetic mathematical constructions. No clinical outcome data are used. Endpoints are computed by the independently audited LP implementation.

### `figures/figure3_simulation.py`

[Open source](figures/figure3_simulation.py) · Supporting module

Figure 3 from verified raw/summary simulation outputs. Shows Monte Carlo uncertainty in certification and interval width. Direct target auditing and one-group projection use different non-target evidence; caption must state that comparison is illustrative, not total-cost matched superiority.

### `figures/figure4_stress_test.py`

[Open source](figures/figure4_stress_test.py) · Supporting module

Version-3 synthetic stress test: assumption failures and interaction. Four heatmap panels show raw-denominator operating characteristics. Only the upper-left design cell within each panel retains the stipulated audit model. This is an illustrative secondary analysis, not a clinical frequency estimate.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--data` | No / positional | default: Path(__file__).resolve().parents[1] / "legacy" / "simulation" /                     "supplement_v4_summary_main_10000.json" |
| `--output` | No / positional | default: Path(__file__).resolve().parent |

### `figures/figureS1_sensitivity_surface.py`

[Open source](figures/figureS1_sensitivity_surface.py) · Supporting module

Supplementary exact-population sensitivity surface for clinical uncertainty. Cells come from the locked V5 rational analytic supplement, not patient data.

### `figures/figureS2_threshold_criterion.py`

[Open source](figures/figureS2_threshold_criterion.py) · Supporting module

Post hoc exact-algebra figure for the fixed-rule threshold corollary. This visualizes stipulated synthetic population bounds and regret, not patient results or a newly selected clinical threshold. Rules and feasible laws remain fixed while the utility threshold changes.

### `legacy/analytic/audit_surface_v5_independent.py`

[Open source](legacy/analytic/audit_surface_v5_independent.py) · Supporting module

Numerical cross-check of V5 exact vertices using a different LP parameterization. This file does not import run_sensitivity_surface_v5.py or reuse its constraints. Its result is a computational cross-check, not external scientific review.

### `legacy/analytic/check_exact.py`

[Open source](legacy/analytic/check_exact.py) · Supporting module

Exact rational audit of Box 1's deterministic identification example. This enumerates vertices of a two-dimensional rational polytope and its intersections with an exact one-stratum audit line. It does not analyse data.

### `legacy/analytic/check_threshold_corollary.py`

[Open source](legacy/analytic/check_threshold_corollary.py) · Executable entry point

Exact check of the post hoc fixed-rule threshold corollary. This reuses only the already established Box 1 sharp endpoints. It does not claim an independent proof of the original feasible-set optimization.

### `legacy/analytic/run_sensitivity_surface_v5.py`

[Open source](legacy/analytic/run_sensitivity_surface_v5.py) · Executable entry point

Exact-rational 56-cell sensitivity surface for Box 1. Design locked in SENSITIVITY_SURFACE_PROTOCOL_V5.md before this code/results. This script only manipulates stipulated synthetic population probabilities.

### `legacy/figures/figure1_design.py`

[Open source](legacy/figures/figure1_design.py) · Supporting module

Figure 1: verifiable, editable vector schematic of the proposed evidence design. All shapes and labels derive from the project's pre-specified analytic/data design. The interval bars are explicitly schematic; no patient results are plotted.

### `legacy/figures/figure2_identification.py`

[Open source](legacy/figures/figure2_identification.py) · Supporting module

Figure 2: exact analytic identified sets for observationally equivalent worlds. The worlds are synthetic mathematical constructions. No clinical outcome data are used. Endpoints are computed by the independently audited LP implementation.

### `legacy/figures/figure3_simulation.py`

[Open source](legacy/figures/figure3_simulation.py) · Supporting module

Figure 3 from verified raw/summary simulation outputs. Shows Monte Carlo uncertainty in certification and interval width. Direct target auditing and one-group projection use different non-target evidence; caption must state that comparison is illustrative, not total-cost matched superiority.

### `legacy/figures/figure4_stress_test.py`

[Open source](legacy/figures/figure4_stress_test.py) · Supporting module

Version-3 synthetic stress test: assumption failures and interaction. Four heatmap panels show raw-denominator operating characteristics. Only the upper-left design cell within each panel retains the stipulated audit model. This is an illustrative secondary analysis, not a clinical frequency estimate.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--data` | No / positional | default: Path(__file__).resolve().parents[1] / "simulation" /                     "supplement_v4_summary_main_10000.json" |
| `--output` | No / positional | default: Path(__file__).resolve().parent |

### `legacy/figures/figureS1_sensitivity_surface.py`

[Open source](legacy/figures/figureS1_sensitivity_surface.py) · Supporting module

Supplementary exact-population sensitivity surface for clinical uncertainty. Cells come from the locked V5 rational analytic supplement, not patient data.

### `legacy/figures/figureS2_threshold_criterion.py`

[Open source](legacy/figures/figureS2_threshold_criterion.py) · Supporting module

Post hoc exact-algebra figure for the fixed-rule threshold corollary. This visualizes stipulated synthetic population bounds and regret, not patient results or a newly selected clinical threshold. Rules and feasible laws remain fixed while the utility threshold changes.

### `legacy/simulation/audit_design_correction_independent.py`

[Open source](legacy/simulation/audit_design_correction_independent.py) · Executable entry point

Independent all-row numerical audit of the post hoc design correction. No production simulation, projection, or correction functions are imported.

### `legacy/simulation/audit_v3_independent.py`

[Open source](legacy/simulation/audit_v3_independent.py) · Supporting module

Independent, standard-library audit of the V3 raw Monte Carlo record. Does not import the simulation or projection implementation. Raises if any record-level flag or aggregated metric conflicts with the saved summary.

### `legacy/simulation/audit_v4_independent.py`

[Open source](legacy/simulation/audit_v4_independent.py) · Executable entry point

Independent row-level V4 audit from saved counts, without simulation imports. Reconstructs seven CP intervals, both 4-variable LPs and all saved results. Also verifies exact V3/V4 row correspondence and immutable V3 hashes.

### `legacy/simulation/run_design_correction_v1.py`

[Open source](legacy/simulation/run_design_correction_v1.py) · Executable entry point

Post hoc mechanism-aware reanalysis of saved V4 stress counts. See DESIGN_CORRECTION_PROTOCOL_V1.md. No new random outcomes are drawn.

### `legacy/simulation/run_study.py`

[Open source](legacy/simulation/run_study.py) · Executable entry point

Internally prespecified finite-sample projection study; see ADEMP_LOCK.md. This is a methodological stress test, not clinical validation or a novel estimator. The retained local run used Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--reps` | Yes |  |
| `--suffix` | Yes |  |

### `legacy/simulation/run_supplement_v2.py`

[Open source](legacy/simulation/run_supplement_v2.py) · Executable entry point

Version-2 post-main-design stress tests and matched target-audit comparator. See SUPPLEMENT_PROTOCOL_V2.md. These results are separate from the frozen main run.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--reps` | Yes |  |
| `--suffix` | Yes |  |

### `legacy/simulation/run_supplement_v3.py`

[Open source](legacy/simulation/run_supplement_v3.py) · Executable entry point

Version-3 2x2 audit-selection/reference-error sensitivity study. See SUPPLEMENT_PROTOCOL_V3.md. All numbers describe stipulated synthetic worlds.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--reps` | Yes |  |
| `--suffix` | Yes |  |

### `legacy/simulation/run_supplement_v4.py`

[Open source](legacy/simulation/run_supplement_v4.py) · Executable entry point

V4 provenance augmentation of the frozen V3 synthetic sensitivity experiment. See SUPPLEMENT_PROTOCOL_V4.md. The original V3 script and outputs are untouched. The added per-row counts permit complete independent recomputation.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--reps` | Yes |  |
| `--suffix` | Yes |  |

### `provenance/october_execution_source/build_figures.py`

[Open source](provenance/october_execution_source/build_figures.py) · Executable entry point

### `provenance/october_execution_source/independent_extension_audit.py`

[Open source](provenance/october_execution_source/independent_extension_audit.py) · Supporting module

Separate reconstruction of E3/E4; imports no author analysis functions.

### `provenance/october_execution_source/run_extensions.py`

[Open source](provenance/october_execution_source/run_extensions.py) · Executable entry point

Revision-stage synthetic analyses. See PROTOCOL.md; no patient data.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--part` | No / positional | default: 'all'; choices: ['joint','budget','all'] |

### `provenance/original_execution_source/run_bounded_correction.py`

[Open source](provenance/original_execution_source/run_bounded_correction.py) · Executable entry point

Frozen synthetic bounded audit correction; see PROTOCOL.md. CPU only.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--reps` | Yes | ; choices: (10,1000) |
| `--label` | Yes | ; choices: ('pilot','main') |

### `revision_20261010/experiments/run_extensions.py`

[Open source](revision_20261010/experiments/run_extensions.py) · Executable entry point

Revision-stage synthetic analyses. See PROTOCOL.md; no patient data.

| Argument | Required | Default / allowed values |
|---|---|---|
| `--part` | No / positional | default: 'all'; choices: ['joint','budget','all'] |

### `revision_20261010/figures/build_figures.py`

[Open source](revision_20261010/figures/build_figures.py) · Executable entry point

### `revision_20261010/verification/independent_extension_audit.py`

[Open source](revision_20261010/verification/independent_extension_audit.py) · Supporting module

Separate reconstruction of E3/E4; imports no author analysis functions.

## Reproduction limits

Absolute paths, frozen-file references and platform-specific assumptions are preserved. Configure a matching local workspace rather than bypassing scientific integrity guards. A repository containing existing scripts is not evidence that every data or environment prerequisite is satisfied. No original clinical results or completed-run artifacts are redistributed.

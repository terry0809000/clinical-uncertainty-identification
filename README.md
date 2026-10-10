# Clinical prediction uncertainty and identification

**Trustworthy AI · Research pipeline · Source-code release**

This project provides the author's existing analysis implementation for clinical prediction uncertainty and identification. Scientific source files are preserved from the local research product. This upload organizes code and documentation; it does not represent new experiments, model retraining or an independently verified reproduction.

## Pipeline

1. **Acquire inputs.** Analytical examples and explicitly synthetic simulations; no patient dataset required for the analytical checks.
2. **Prepare the environment.** Use the project-specific dependency files listed below in an isolated environment. Preserve the versions and seeds recorded in the source release.
3. **Prepare data and run the model or analysis.** Follow the original execution notes and entry points below, preserving time order, participant grouping and fixed evaluation boundaries.
4. **Evaluate and summarize.** Use the included evaluation/analysis scripts. Write generated data, predictions, fitted objects and caches to a private working copy.

This is a source-only repository. Raw data, individual-level derivatives and fitted checkpoints are not bundled. Some saved-evidence reproduction commands in the original notes require artifacts from the original source release; obtain those independently where authorized.

## Environment

- [`legacy/requirements.txt`](legacy/requirements.txt)
- [`requirements_extensions.txt`](requirements_extensions.txt)
- [`requirements_revision.txt`](requirements_revision.txt)
- [`revision_20261010/experiments/environment.json`](revision_20261010/experiments/environment.json)

## Execution

See the [pipeline and command reference](PIPELINE.md) for entry points, required arguments and the source-specific input contract.

Run commands from the repository root unless the original notes specify another working directory. Use a fresh copy and configure local data paths before running.

### Original execution notes

- [README.md](docs/original/README.md)
- [README.md](docs/original/extensions/bounded_correction/README.md)
- [README.txt](docs/original/legacy/README.txt)

### Scientific source inventory

<details>
<summary>Browse the preserved scientific implementation</summary>

- [`extensions/bounded_correction/audit_bounded_correction.py`](extensions/bounded_correction/audit_bounded_correction.py)
- [`extensions/bounded_correction/run_bounded_correction.py`](extensions/bounded_correction/run_bounded_correction.py)
- [`extensions/group_channels/group_channels.py`](extensions/group_channels/group_channels.py)
- [`figures/build_figure5_and_exports.py`](figures/build_figure5_and_exports.py)
- [`figures/figure1_evidence_design.py`](figures/figure1_evidence_design.py)
- [`figures/figure2_identification.py`](figures/figure2_identification.py)
- [`figures/figure3_simulation.py`](figures/figure3_simulation.py)
- [`figures/figure4_stress_test.py`](figures/figure4_stress_test.py)
- [`figures/figureS1_sensitivity_surface.py`](figures/figureS1_sensitivity_surface.py)
- [`figures/figureS2_threshold_criterion.py`](figures/figureS2_threshold_criterion.py)
- [`legacy/analytic/audit_surface_v5_independent.py`](legacy/analytic/audit_surface_v5_independent.py)
- [`legacy/analytic/check_exact.py`](legacy/analytic/check_exact.py)
- [`legacy/analytic/check_threshold_corollary.py`](legacy/analytic/check_threshold_corollary.py)
- [`legacy/analytic/run_sensitivity_surface_v5.py`](legacy/analytic/run_sensitivity_surface_v5.py)
- [`legacy/figures/figure1_design.py`](legacy/figures/figure1_design.py)
- [`legacy/figures/figure2_identification.py`](legacy/figures/figure2_identification.py)
- [`legacy/figures/figure3_simulation.py`](legacy/figures/figure3_simulation.py)
- [`legacy/figures/figure4_stress_test.py`](legacy/figures/figure4_stress_test.py)
- [`legacy/figures/figureS1_sensitivity_surface.py`](legacy/figures/figureS1_sensitivity_surface.py)
- [`legacy/figures/figureS2_threshold_criterion.py`](legacy/figures/figureS2_threshold_criterion.py)
- [`legacy/simulation/audit_design_correction_independent.py`](legacy/simulation/audit_design_correction_independent.py)
- [`legacy/simulation/audit_v3_independent.py`](legacy/simulation/audit_v3_independent.py)
- [`legacy/simulation/audit_v4_independent.py`](legacy/simulation/audit_v4_independent.py)
- [`legacy/simulation/run_design_correction_v1.py`](legacy/simulation/run_design_correction_v1.py)
- [`legacy/simulation/run_study.py`](legacy/simulation/run_study.py)
- [`legacy/simulation/run_supplement_v2.py`](legacy/simulation/run_supplement_v2.py)
- [`legacy/simulation/run_supplement_v3.py`](legacy/simulation/run_supplement_v3.py)
- [`legacy/simulation/run_supplement_v4.py`](legacy/simulation/run_supplement_v4.py)
- [`provenance/october_execution_source/build_figures.py`](provenance/october_execution_source/build_figures.py)
- [`provenance/october_execution_source/independent_extension_audit.py`](provenance/october_execution_source/independent_extension_audit.py)
- [`provenance/october_execution_source/run_extensions.py`](provenance/october_execution_source/run_extensions.py)
- [`provenance/original_execution_source/run_bounded_correction.py`](provenance/original_execution_source/run_bounded_correction.py)
- [`revision_20261010/experiments/run_extensions.py`](revision_20261010/experiments/run_extensions.py)
- [`revision_20261010/figures/build_figures.py`](revision_20261010/figures/build_figures.py)
- [`revision_20261010/verification/independent_extension_audit.py`](revision_20261010/verification/independent_extension_audit.py)

</details>

## Repository contents and rights

Only scientific code, dependencies, configuration and execution documentation are included. Manuscripts, submission documents, raw participant records, trained models, credentials, caches and software environments are excluded. Notebook cell outputs and execution counters are cleared for publication; code cells are retained.

Existing source and third-party license notices govern their respective files. No blanket license is invented for material whose original release did not specify one. Original dataset licenses and access requirements continue to apply.

Research software is not a clinical deployment claim. Negative findings and the source release's interpretation limits remain applicable. No tests or experiments were run during this upload, as requested by the author.

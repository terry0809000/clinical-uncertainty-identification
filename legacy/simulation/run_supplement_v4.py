"""V4 provenance augmentation of the frozen V3 synthetic sensitivity experiment.

See SUPPLEMENT_PROTOCOL_V4.md. The original V3 script and outputs are untouched.
The added per-row counts permit complete independent recomputation.
"""

from __future__ import annotations

import argparse
import csv
import json
from math import sqrt
from pathlib import Path

import numpy as np

from run_study import (
    ALPHA, N_SOURCE, N_TARGET, SOURCE_P, SOURCE_SE, SOURCE_SP,
    World, cp_interval, delta, envelope, exact_population_interval, make_worlds,
)
from run_supplement_v2 import new_worlds, upper_binomial_95
from run_supplement_v3 import (
    SEED, N_AUDIT, LAMBDAS, REFERENCE_CHANNELS, audit_positive_probability,
)

HERE = Path(__file__).resolve().parent
PRIMITIVES = (
    "source_gold_plus", "source_gold_minus", "source_se", "source_sp",
    "target_proxy_plus", "target_proxy_minus", "target_audit_recorded",
)
BASE_FIELDS = (
    "world", "lambda_proxy_positive", "ref_se", "ref_sp", "rep",
    "audit_selected_true_risk", "audit_recorded_positive_rate",
    "lower", "upper", "true_delta", "true_covered",
    "correct_model_set_covered", "certification", "correct_sign",
    "wrong_sign", "indeterminate", "infeasible", "width",
)
ADDED_FIELDS = tuple(
    f"{name}_{field}"
    for name in PRIMITIVES
    for field in ("success", "denominator", "lower", "upper")
)


def sample_logged(rng: np.random.Generator, world: World, recorded_rate: float):
    """Match V3 RNG draw order; save all seven primitive sufficient statistics."""
    k = [int(rng.binomial(N_SOURCE, p)) for p in SOURCE_P]
    se_yes = sum(int(rng.binomial(ki, SOURCE_SE)) for ki in k)
    sp_yes = sum(int(rng.binomial(N_SOURCE - ki, SOURCE_SP)) for ki in k)
    q_plus = int(rng.binomial(N_TARGET, world.q_plus))
    q_minus = int(rng.binomial(N_TARGET, world.q_minus))
    # Reproduce the discarded draw in V3 sample_primitives, without using it.
    _unused_unbiased_audit = int(rng.binomial(N_AUDIT, world.p_plus))
    recorded = int(rng.binomial(N_AUDIT, recorded_rate))
    observations = (
        (k[0], N_SOURCE),
        (k[1], N_SOURCE),
        (se_yes, sum(k)),
        (sp_yes, 2 * N_SOURCE - sum(k)),
        (q_plus, N_TARGET),
        (q_minus, N_TARGET),
        (recorded, N_AUDIT),
    )
    intervals = tuple(cp_interval(x, n, ALPHA / 7) for x, n in observations)
    evidence = {
        f"{name}_{part}": value
        for name, (x, n), (lo, hi) in zip(PRIMITIVES, observations, intervals)
        for part, value in (
            ("success", x), ("denominator", n), ("lower", lo), ("upper", hi)
        )
    }
    return intervals, evidence


def run(reps: int, suffix: str):
    _, b, _ = make_worlds()
    c, _ = new_worlds()
    worlds = (b, c)
    designs = [
        (world, lam, ref_se, ref_sp)
        for world in worlds
        for lam in LAMBDAS
        for ref_se, ref_sp in REFERENCE_CHANNELS
    ]
    assert len(designs) == 8
    seeds = np.random.SeedSequence(SEED).spawn(len(designs))
    raw_path = HERE / f"supplement_v4_raw_{suffix}.csv"
    summary_path = HERE / f"supplement_v4_summary_{suffix}.json"
    summary = {
        "protocol": "SUPPLEMENT_PROTOCOL_V4.md",
        "seed": SEED,
        "repetitions_per_cell": reps,
        "n_target_reference_audit": N_AUDIT,
        "cells": [],
        "interactions": [],
    }
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_FIELDS + ADDED_FIELDS)
        writer.writeheader()
        for cell_index, (world, lam, ref_se, ref_sp) in enumerate(designs):
            rng = np.random.default_rng(seeds[cell_index])
            selected_risk, recorded_rate = audit_positive_probability(
                world, lam, ref_se, ref_sp
            )
            true_delta = delta(world.p_plus, world.p_minus)
            pop_interval = exact_population_interval(world, N_AUDIT)
            valid = lam == 1.0 and ref_se == 1.0 and ref_sp == 1.0
            counts = {
                name: 0 for name in (
                    "true_covered", "correct_model_set_covered",
                    "correct_sign", "wrong_sign", "indeterminate", "infeasible",
                )
            }
            widths = []
            for rep in range(reps):
                intervals, evidence = sample_logged(rng, world, recorded_rate)
                result = envelope(
                    intervals[:2], intervals[2], intervals[3],
                    intervals[4:6], intervals[6],
                )
                infeasible = result is None
                if infeasible:
                    lower = upper = width = float("nan")
                    certification = "infeasible"
                else:
                    lower, upper = result
                    width = upper - lower
                    widths.append(width)
                    certification = (
                        "positive" if lower > 1e-12
                        else "negative" if upper < -1e-12
                        else "indeterminate"
                    )
                covered = (
                    not infeasible
                    and lower - 1e-9 <= true_delta <= upper + 1e-9
                )
                set_covered = (
                    not infeasible
                    and lower - 1e-9 <= pop_interval[0]
                    and upper + 1e-9 >= pop_interval[1]
                )
                true_sign = (
                    "positive" if true_delta > 1e-12
                    else "negative" if true_delta < -1e-12
                    else "tie"
                )
                correct = certification == true_sign
                wrong = (
                    certification in ("positive", "negative")
                    and certification != true_sign
                )
                flags = {
                    "true_covered": covered,
                    "correct_model_set_covered": set_covered,
                    "correct_sign": correct,
                    "wrong_sign": wrong,
                    "indeterminate": certification == "indeterminate",
                    "infeasible": infeasible,
                }
                for name, flag in flags.items():
                    counts[name] += int(flag)
                writer.writerow({
                    "world": world.name,
                    "lambda_proxy_positive": lam,
                    "ref_se": ref_se,
                    "ref_sp": ref_sp,
                    "rep": rep,
                    "audit_selected_true_risk": selected_risk,
                    "audit_recorded_positive_rate": recorded_rate,
                    "lower": lower,
                    "upper": upper,
                    "true_delta": true_delta,
                    "true_covered": int(covered),
                    "correct_model_set_covered": int(set_covered),
                    "certification": certification,
                    "correct_sign": int(correct),
                    "wrong_sign": int(wrong),
                    "indeterminate": int(certification == "indeterminate"),
                    "infeasible": int(infeasible),
                    "width": width,
                    **evidence,
                })
            metrics = {
                name: {
                    "count": count,
                    "proportion": count / reps,
                    "mcse": sqrt((count / reps) * (1 - count / reps) / reps),
                    "one_sided_95_upper": upper_binomial_95(count, reps),
                }
                for name, count in counts.items()
            }
            metrics["mean_width_when_feasible"] = (
                float(np.mean(widths)) if widths else None
            )
            summary["cells"].append({
                "world": world.name,
                "lambda_proxy_positive": lam,
                "reference_se": ref_se,
                "reference_sp": ref_sp,
                "assumptions_valid": valid,
                "selected_true_risk": selected_risk,
                "recorded_audit_positive_rate": recorded_rate,
                "true_delta": true_delta,
                "correct_model_population_set": pop_interval,
                "metrics": metrics,
            })
            summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
            print(
                world.name, "lambda", lam, "reference", ref_se, ref_sp,
                "coverage", metrics["true_covered"]["proportion"],
                "wrong", metrics["wrong_sign"]["proportion"],
                "infeasible", metrics["infeasible"]["proportion"],
                flush=True,
            )
        for world in worlds:
            cells = [cell for cell in summary["cells"] if cell["world"] == world.name]
            for outcome in ("true_covered", "wrong_sign"):
                rates, variances = {}, {}
                for cell in cells:
                    key = (
                        int(cell["lambda_proxy_positive"] == 4.0),
                        int(cell["reference_se"] < 1.0),
                    )
                    p = cell["metrics"][outcome]["proportion"]
                    rates[key] = p
                    variances[key] = p * (1 - p) / reps
                interaction = (
                    rates[1, 1] - rates[1, 0] - rates[0, 1] + rates[0, 0]
                )
                mcse = sqrt(sum(variances.values()))
                summary["interactions"].append({
                    "world": world.name,
                    "outcome": outcome,
                    "probability_scale_interaction": interaction,
                    "mcse": mcse,
                    "descriptive_95_monte_carlo_interval": [
                        interaction - 1.96 * mcse,
                        interaction + 1.96 * mcse,
                    ],
                })
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("wrote", raw_path, summary_path, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument("--suffix", required=True)
    args = parser.parse_args()
    if args.reps <= 0 or not args.suffix.replace("_", "").isalnum():
        parser.error("Positive repetitions and alphanumeric/underscore suffix required")
    run(args.reps, args.suffix)

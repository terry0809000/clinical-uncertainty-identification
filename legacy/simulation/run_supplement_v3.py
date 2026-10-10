"""Version-3 2x2 audit-selection/reference-error sensitivity study.

See SUPPLEMENT_PROTOCOL_V3.md. All numbers describe stipulated synthetic worlds.
"""

from __future__ import annotations

import argparse
import csv
import json
from math import sqrt
from pathlib import Path

import numpy as np

from run_study import (
    ALPHA,
    World,
    cp_interval,
    delta,
    envelope,
    exact_population_interval,
    make_worlds,
    sample_primitives,
)
from run_supplement_v2 import new_worlds, upper_binomial_95

HERE = Path(__file__).resolve().parent
SEED = 20260926
N_AUDIT = 800
LAMBDAS = (1.0, 4.0)
REFERENCE_CHANNELS = ((1.0, 1.0), (0.80, 0.95))


def audit_positive_probability(world: World, lam: float, ref_se: float, ref_sp: float):
    """Proxy-weighted true risk, then imperfect-reference positive probability."""
    p, se, sp = world.p_plus, world.se_plus, world.sp_plus
    p11 = p * se
    p10 = p * (1 - se)
    p01 = (1 - p) * (1 - sp)
    p00 = (1 - p) * sp
    q = p11 + p01
    assert abs(q - world.q_plus) < 1e-12
    normalizer = lam * (p11 + p01) + p10 + p00
    selected_true_risk = (lam * p11 + p10) / normalizer
    reference_positive_rate = (
        ref_se * selected_true_risk + (1 - ref_sp) * (1 - selected_true_risk)
    )
    assert 0 <= reference_positive_rate <= 1
    return selected_true_risk, reference_positive_rate


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
    raw_path = HERE / f"supplement_v3_raw_{suffix}.csv"
    summary_path = HERE / f"supplement_v3_summary_{suffix}.json"
    seeds = np.random.SeedSequence(SEED).spawn(len(designs))
    summary = {
        "protocol": "SUPPLEMENT_PROTOCOL_V3.md",
        "seed": SEED,
        "repetitions_per_cell": reps,
        "n_target_reference_audit": N_AUDIT,
        "cells": [],
        "interactions": [],
    }
    fieldnames = (
        "world", "lambda_proxy_positive", "ref_se", "ref_sp", "rep",
        "audit_selected_true_risk", "audit_recorded_positive_rate",
        "lower", "upper", "true_delta", "true_covered",
        "correct_model_set_covered", "certification", "correct_sign",
        "wrong_sign", "indeterminate", "infeasible", "width",
    )
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
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
                name: 0
                for name in (
                    "true_covered", "correct_model_set_covered",
                    "correct_sign", "wrong_sign", "indeterminate", "infeasible",
                )
            }
            widths = []
            for rep in range(reps):
                source_p, source_se, source_sp, q, _ = sample_primitives(
                    rng, world, N_AUDIT
                )
                positive = int(rng.binomial(N_AUDIT, recorded_rate))
                audit = cp_interval(positive, N_AUDIT, ALPHA / 7)
                interval = envelope(source_p, source_se, source_sp, q, audit)
                infeasible = interval is None
                if infeasible:
                    lower = upper = float("nan")
                    certification = "infeasible"
                    width = float("nan")
                else:
                    lower, upper = interval
                    width = upper - lower
                    widths.append(width)
                    certification = (
                        "positive" if lower > 1e-12
                        else "negative" if upper < -1e-12
                        else "indeterminate"
                    )
                covered = not infeasible and lower - 1e-9 <= true_delta <= upper + 1e-9
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
                wrong = certification in ("positive", "negative") and certification != true_sign
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
                writer.writerow(
                    {
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
                    }
                )
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
            summary["cells"].append(
                {
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
                }
            )
            print(
                world.name, "lambda", lam, "reference", ref_se, ref_sp,
                "coverage", metrics["true_covered"]["proportion"],
                "wrong", metrics["wrong_sign"]["proportion"],
                "infeasible", metrics["infeasible"]["proportion"],
                flush=True,
            )
            summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        for world in worlds:
            cells = [
                cell for cell in summary["cells"] if cell["world"] == world.name
            ]
            for outcome in ("true_covered", "wrong_sign"):
                rates = {}
                variances = {}
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
                summary["interactions"].append(
                    {
                        "world": world.name,
                        "outcome": outcome,
                        "probability_scale_interaction": interaction,
                        "mcse": mcse,
                        "descriptive_95_monte_carlo_interval": [
                            interaction - 1.96 * mcse,
                            interaction + 1.96 * mcse,
                        ],
                    }
                )
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("wrote", raw_path, summary_path, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument("--suffix", required=True)
    arguments = parser.parse_args()
    if arguments.reps <= 0 or not arguments.suffix.replace("_", "").isalnum():
        parser.error("Positive repetitions and alphanumeric/underscore suffix required")
    run(arguments.reps, arguments.suffix)

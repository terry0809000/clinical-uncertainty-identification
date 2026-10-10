"""Post hoc mechanism-aware reanalysis of saved V4 stress counts.

See DESIGN_CORRECTION_PROTOCOL_V1.md. No new random outcomes are drawn.
"""

from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from scipy.stats import beta

from run_study import envelope

HERE = Path(__file__).resolve().parent
INPUT = HERE / "supplement_v4_raw_main_10000.csv"
OUTPUT = HERE / "design_correction_v1_raw.csv"
SUMMARY = HERE / "design_correction_v1_summary.json"
PROTOCOL = HERE / "DESIGN_CORRECTION_PROTOCOL_V1.md"
ALPHA = 0.05 / 7
N = 10000
WORLD_TRUE_P = {"B_negative": 19 / 60, "C_interior_tie": 7 / 20}
WORLD_SET_CONTROL = {
    "B_negative": (-1 / 72, -1 / 72),
    "C_interior_tie": (-1 / 152, 1 / 328),
}
PRIMITIVES = (
    "source_gold_plus", "source_gold_minus", "source_se", "source_sp",
    "target_proxy_plus", "target_proxy_minus", "target_audit_recorded",
)
FIELDS = (
    "world", "lambda_proxy_positive", "ref_se", "ref_sp", "rep",
    "corrected_p_lower", "corrected_p_upper", "lower", "upper", "true_delta",
    "true_covered", "correct_sign", "wrong_sign", "indeterminate",
    "infeasible", "width", "control_set_covered",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def cp(x: int, n: int) -> tuple[float, float]:
    if n == 0:
        return 0.0, 1.0
    lo = 0.0 if x == 0 else float(beta.ppf(ALPHA / 2, x, n - x + 1))
    hi = 1.0 if x == n else float(beta.ppf(1 - ALPHA / 2, x + 1, n - x))
    return lo, hi


def mc_interval(events: int, trials: int) -> tuple[float, float]:
    lo = 0.0 if events == 0 else float(beta.ppf(0.025, events, trials - events + 1))
    hi = 1.0 if events == trials else float(beta.ppf(0.975, events + 1, trials - events))
    return lo, hi


def corrected_p_interval(
    audit: tuple[float, float],
    q: tuple[float, float],
    target_se: tuple[float, float],
    lam: float,
    ref_se: float,
    ref_sp: float,
) -> tuple[float, float] | None:
    assert lam >= 1 and ref_se + ref_sp > 1
    a_low = max(audit[0], 1 - ref_sp)
    a_high = min(audit[1], ref_se)
    if a_low > a_high:
        return None
    gain = ref_se + ref_sp - 1
    r_low = (a_low + ref_sp - 1) / gain
    r_high = (a_high + ref_sp - 1) / gain
    k = lam - 1
    lo = r_low * (1 + k * q[0]) / (1 + k * target_se[1])
    hi = r_high * (1 + k * q[1]) / (1 + k * target_se[0])
    lo, hi = max(0.0, lo), min(1.0, hi)
    if lo > hi + 1e-12:
        return None
    return lo, hi


def main() -> None:
    counts = defaultdict(lambda: defaultdict(int))
    widths = defaultdict(list)
    max_cp_discrepancy = 0.0
    max_control_discrepancy = 0.0
    rows = 0
    with INPUT.open(newline="", encoding="utf-8") as fin, OUTPUT.open(
        "w", newline="", encoding="utf-8"
    ) as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=FIELDS)
        writer.writeheader()
        for row in reader:
            world = row["world"]
            lam = float(row["lambda_proxy_positive"])
            ref_se = float(row["ref_se"])
            ref_sp = float(row["ref_sp"])
            key = (world, lam, ref_se, ref_sp)
            assert world in WORLD_TRUE_P
            intervals = {}
            for primitive in PRIMITIVES:
                this = cp(
                    int(row[f"{primitive}_success"]),
                    int(row[f"{primitive}_denominator"]),
                )
                saved = (
                    float(row[f"{primitive}_lower"]),
                    float(row[f"{primitive}_upper"]),
                )
                max_cp_discrepancy = max(
                    max_cp_discrepancy,
                    abs(this[0] - saved[0]), abs(this[1] - saved[1]),
                )
                intervals[primitive] = this
            source_se = intervals["source_se"]
            target_se = (max(0.0, source_se[0] - 0.05),
                         min(1.0, source_se[1] + 0.05))
            p_bound = corrected_p_interval(
                intervals["target_audit_recorded"],
                intervals["target_proxy_plus"], target_se,
                lam, ref_se, ref_sp,
            )
            projected = None if p_bound is None else envelope(
                (intervals["source_gold_plus"], intervals["source_gold_minus"]),
                intervals["source_se"], intervals["source_sp"],
                (intervals["target_proxy_plus"], intervals["target_proxy_minus"]),
                p_bound,
            )
            true_delta = float(row["true_delta"])
            if projected is None:
                low = high = width = None
                covered = correct = wrong = indeterminate = set_covered = 0
                infeasible = 1
            else:
                low, high = projected
                assert low <= high + 1e-10
                width = high - low
                covered = int(low - 1e-10 <= true_delta <= high + 1e-10)
                wrong = int((low > 1e-10 and true_delta <= 0) or
                            (high < -1e-10 and true_delta >= 0))
                correct = int((low > 1e-10 and true_delta > 0) or
                              (high < -1e-10 and true_delta < 0))
                indeterminate = int(not (wrong or correct))
                infeasible = 0
                widths[key].append(width)
                set_covered = 0
                if lam == 1 and ref_se == 1 and ref_sp == 1:
                    target_set = WORLD_SET_CONTROL[world]
                    set_covered = int(low - 1e-10 <= target_set[0] and
                                      target_set[1] <= high + 1e-10)
                    max_control_discrepancy = max(
                        max_control_discrepancy,
                        abs(low - float(row["lower"])),
                        abs(high - float(row["upper"])),
                    )
            counts[key]["n"] += 1
            for metric, value in (
                ("true_covered", covered), ("correct_sign", correct),
                ("wrong_sign", wrong), ("indeterminate", indeterminate),
                ("infeasible", infeasible),
                ("control_set_covered", set_covered),
                ("naive_true_covered", int(row["true_covered"])),
                ("naive_wrong_sign", int(row["wrong_sign"])),
            ):
                counts[key][metric] += value
            writer.writerow({
                "world": world, "lambda_proxy_positive": lam,
                "ref_se": ref_se, "ref_sp": ref_sp, "rep": row["rep"],
                "corrected_p_lower": "" if p_bound is None else p_bound[0],
                "corrected_p_upper": "" if p_bound is None else p_bound[1],
                "lower": "" if low is None else low,
                "upper": "" if high is None else high,
                "true_delta": true_delta,
                "true_covered": covered, "correct_sign": correct,
                "wrong_sign": wrong, "indeterminate": indeterminate,
                "infeasible": infeasible,
                "width": "" if width is None else width,
                "control_set_covered": set_covered if lam == 1 and ref_se == 1 and ref_sp == 1 else "",
            })
            rows += 1
            if rows % 10000 == 0:
                print(f"processed {rows}", flush=True)
    assert rows == 80000 and len(counts) == 8
    assert max_cp_discrepancy < 1e-12, max_cp_discrepancy
    assert max_control_discrepancy < 1e-10, max_control_discrepancy
    cells = []
    for key in sorted(counts):
        c = counts[key]
        assert c["n"] == N
        assert c["true_covered"] + c["infeasible"] <= N
        assert c["correct_sign"] + c["wrong_sign"] + c["indeterminate"] + c["infeasible"] == N
        metrics = {}
        for metric in ("true_covered", "correct_sign", "wrong_sign", "indeterminate",
                       "infeasible", "naive_true_covered", "naive_wrong_sign"):
            value = c[metric]
            low, high = mc_interval(value, N)
            metrics[metric] = {"count": value, "proportion": value / N,
                               "mc_95_interval": [low, high]}
        if key[1:] == (1.0, 1.0, 1.0):
            value = c["control_set_covered"]
            metrics["control_set_covered"] = {
                "count": value, "proportion": value / N,
                "mc_95_interval": list(mc_interval(value, N)),
            }
        metrics["mean_width_feasible_only"] = (
            sum(widths[key]) / len(widths[key]) if widths[key] else None
        )
        cells.append({"world": key[0], "lambda_proxy_positive": key[1],
                      "ref_se": key[2], "ref_sp": key[3], "metrics": metrics})
        print(key, "coverage", c["true_covered"] / N,
              "wrong_sign", c["wrong_sign"] / N,
              "infeasible", c["infeasible"] / N, flush=True)
    result = {
        "status": "POST_HOC_EXPLORATORY_MECHANISM_KNOWN",
        "protocol": PROTOCOL.name, "protocol_sha256": sha(PROTOCOL),
        "input": INPUT.name, "input_sha256": sha(INPUT),
        "output": OUTPUT.name, "output_sha256": sha(OUTPUT),
        "script_sha256": sha(Path(__file__)),
        "n_total": rows, "n_per_cell": N,
        "max_cp_saved_endpoint_discrepancy": max_cp_discrepancy,
        "max_control_naive_endpoint_discrepancy": max_control_discrepancy,
        "cells": cells,
        "limits": "No new simulated patients; known selection and recording channels; mechanism-specific sharp population sets not derived outside control; same-family analysis pending independent audit.",
    }
    SUMMARY.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Version-2 post-main-design stress tests and matched target-audit comparator.

See SUPPLEMENT_PROTOCOL_V2.md. These results are separate from the frozen main run.
"""

from __future__ import annotations

import argparse
import csv
import json
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.stats import beta

from run_study import (ALPHA, AUDIT_SIZES, SOURCE_P, SOURCE_SE, SOURCE_SP,
                       DRIFT_P, DRIFT_SE, DRIFT_SP, ETA, World,
                       cp_interval, delta, envelope, exact_population_interval,
                       make_worlds, sample_primitives)

HERE = Path(__file__).resolve().parent
SEED_V2 = 20260925


def new_worlds():
    c = World("C_interior_tie", 0.35, 0.50, 0.31, 0.40,
              0.70, 0.90, 0.70, 0.90, True)
    d = World("D_heterogeneous_only", 0.46, 0.49, 0.3114, 0.3536,
              0.63, 0.96, 0.68, 0.96, False)
    for w in (c, d):
        for p, q, se, sp, p0 in ((w.p_plus, w.q_plus, w.se_plus, w.sp_plus, SOURCE_P[0]),
                                 (w.p_minus, w.q_minus, w.se_minus, w.sp_minus, SOURCE_P[1])):
            assert abs(q - (se * p + (1 - sp) * (1 - p))) < 1e-12
            assert abs(p - p0) <= DRIFT_P + 1e-12
            assert abs(se - SOURCE_SE) <= DRIFT_SE + 1e-12
            assert abs(sp - SOURCE_SP) <= DRIFT_SP + 1e-12
            assert se + sp - 1 >= ETA - 1e-12
    assert abs(delta(c.p_plus, c.p_minus)) < 1e-12
    assert abs(delta(d.p_plus, d.p_minus) - 23 / 240) < 1e-12
    return c, d


def upper_binomial_95(events, reps):
    return 1.0 if events == reps else float(beta.ppf(0.95, events + 1, reps - events))


def direct_audit_interval(rng, w, n):
    np_, nm = n // 2, n - n // 2
    kp = int(rng.binomial(np_, w.p_plus))
    km = int(rng.binomial(nm, w.p_minus))
    # Two primitive intervals, Bonferroni alpha/2 each.
    lp, up = cp_interval(kp, np_, ALPHA / 2)
    lm, um = cp_interval(km, nm, ALPHA / 2)
    return delta(lp, um), delta(up, lm)


def summarise(flags, widths, n):
    ans = {key: {"count": count, "proportion": count / n,
                 "mcse": sqrt((count / n) * (1 - count / n) / n),
                 "one_sided_95_upper": upper_binomial_95(count, n)}
           for key, count in flags.items()}
    ans["mean_interval_width"] = float(np.mean(widths)) if widths else None
    ans["median_interval_width"] = float(np.median(widths)) if widths else None
    return ans


def run(reps: int, suffix: str):
    c, d = new_worlds()
    a, b, _ = make_worlds()
    world_by_name = {w.name: w for w in (a, b, c, d)}
    seeds = iter(np.random.SeedSequence(SEED_V2).spawn(20))
    raw_path = HERE / f"supplement_v2_raw_{suffix}.csv"
    summary_path = HERE / f"supplement_v2_summary_{suffix}.json"
    summary = {"protocol": "SUPPLEMENT_PROTOCOL_V2.md", "seed": SEED_V2,
               "reps_per_cell": reps, "cells": []}
    with raw_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=("method", "world", "n_audit", "rep", "lower", "upper",
                                               "true_delta", "true_covered", "set_covered", "certification",
                                               "correct_sign", "wrong_sign", "infeasible", "width"))
        writer.writeheader()
        designs = [("one_group_projection", w, n) for w in (c, d) for n in AUDIT_SIZES]
        designs += [("direct_two_group_audit", w, n) for w in (a, b, c, d) for n in AUDIT_SIZES[1:]]
        assert len(designs) == 20
        for method, w, n in designs:
            rng = np.random.default_rng(next(seeds))
            truth = delta(w.p_plus, w.p_minus)
            population_set = exact_population_interval(w, n) if method == "one_group_projection" else None
            flags = {key: 0 for key in ("true_covered", "set_covered", "correct_sign", "wrong_sign",
                                        "infeasible", "indeterminate")}
            widths = []
            for rep in range(reps):
                if method == "one_group_projection":
                    interval = envelope(*sample_primitives(rng, w, n))
                else:
                    interval = direct_audit_interval(rng, w, n)
                infeasible = interval is None
                if infeasible:
                    low = high = float("nan")
                    cert = "infeasible"
                else:
                    low, high = interval
                    widths.append(high - low)
                    cert = "positive" if low > 1e-12 else "negative" if high < -1e-12 else "indeterminate"
                true_covered = not infeasible and low - 1e-9 <= truth <= high + 1e-9
                set_covered = (None if method != "one_group_projection" or not w.assumed_valid or
                               population_set is None else
                               (not infeasible and low - 1e-9 <= population_set[0] and
                                high + 1e-9 >= population_set[1]))
                true_sign = "positive" if truth > 1e-12 else "negative" if truth < -1e-12 else "tie"
                correct = cert == true_sign
                wrong = cert in ("positive", "negative") and cert != true_sign
                values = {"true_covered": true_covered, "set_covered": set_covered,
                          "correct_sign": correct, "wrong_sign": wrong, "infeasible": infeasible,
                          "indeterminate": cert == "indeterminate"}
                for key, val in values.items():
                    if val is not None:
                        flags[key] += int(val)
                writer.writerow({"method": method, "world": w.name, "n_audit": n, "rep": rep,
                                 "lower": low, "upper": high, "true_delta": truth,
                                 "true_covered": int(true_covered),
                                 "set_covered": "" if set_covered is None else int(set_covered),
                                 "certification": cert, "correct_sign": int(correct), "wrong_sign": int(wrong),
                                 "infeasible": int(infeasible), "width": high - low if not infeasible else float("nan")})
            metrics = summarise(flags, widths, reps)
            if method != "one_group_projection" or not w.assumed_valid or population_set is None:
                metrics["set_covered"] = None
            summary["cells"].append({"method": method, "world": w.name,
                                     "assumptions_valid": w.assumed_valid if method == "one_group_projection" else True,
                                     "n_audit": n, "true_delta": truth,
                                     "population_set": population_set, "metrics": metrics})
            print(method, w.name, n, "coverage", metrics["true_covered"]["proportion"],
                  "wrong", metrics["wrong_sign"]["proportion"],
                  "infeasible", metrics["infeasible"]["proportion"], flush=True)
            summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("wrote", raw_path, summary_path, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument("--suffix", required=True)
    args = parser.parse_args()
    if args.reps <= 0 or not args.suffix.replace("_", "").isalnum():
        parser.error("Positive repetitions and alphanumeric suffix required")
    run(args.reps, args.suffix)

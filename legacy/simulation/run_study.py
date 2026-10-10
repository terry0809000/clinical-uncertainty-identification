"""Internally prespecified finite-sample projection study; see ADEMP_LOCK.md.

This is a methodological stress test, not clinical validation or a novel estimator.
The retained local run used Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.stats import beta


HERE = Path(__file__).resolve().parent
SEED = 20260924
ALPHA = 0.05
T = 0.20
WEIGHTS = (2 / 3, 1 / 3)
SOURCE_P = (0.40, 0.55)
SOURCE_SE, SOURCE_SP = 0.67, 0.93
DRIFT_P, DRIFT_SE, DRIFT_SP = 0.10, 0.05, 0.05
ETA = 0.50
N_SOURCE = 1000
N_TARGET = 1000
AUDIT_SIZES = (0, 50, 200, 800)


@dataclass(frozen=True)
class World:
    name: str
    p_plus: float
    p_minus: float
    q_plus: float
    q_minus: float
    se_plus: float
    sp_plus: float
    se_minus: float
    sp_minus: float
    assumed_valid: bool


def make_worlds() -> tuple[World, ...]:
    a_p, a_m = 29 / 60, 19 / 30
    b_p, b_m = 19 / 60, 7 / 15
    v_se_p = 0.90
    v_sp_p = 1 - (0.31 - v_se_p * b_p) / (1 - b_p)
    return (
        World("A_positive", a_p, a_m, 0.31, 0.40, 0.62, 0.98, 0.62, 0.98, True),
        World("B_negative", b_p, b_m, 0.31, 0.40, 0.72, 0.88, 0.72, 0.88, True),
        World("V_differential_error", b_p, b_m, 0.31, 0.30,
              v_se_p, v_sp_p, 0.62, 0.98, False),
    )


def delta(p_plus: float, p_minus: float) -> float:
    wp, wm = WEIGHTS
    return (wp * p_plus - wm * p_minus - T * (wp - wm)) / (1 - T)


def validate_worlds(worlds: tuple[World, ...]) -> None:
    for w in worlds:
        for x in (w.p_plus, w.p_minus, w.q_plus, w.q_minus,
                  w.se_plus, w.sp_plus, w.se_minus, w.sp_minus):
            assert 0 <= x <= 1, (w.name, x)
        for p, q, se, sp in ((w.p_plus, w.q_plus, w.se_plus, w.sp_plus),
                             (w.p_minus, w.q_minus, w.se_minus, w.sp_minus)):
            assert abs(q - (se * p + (1 - sp) * (1 - p))) < 1e-12, w.name
        if w.assumed_valid:
            for p, p0 in zip((w.p_plus, w.p_minus), SOURCE_P):
                assert abs(p - p0) <= DRIFT_P + 1e-12
            for se, sp in ((w.se_plus, w.sp_plus), (w.se_minus, w.sp_minus)):
                assert abs(se - SOURCE_SE) <= DRIFT_SE + 1e-12
                assert abs(sp - SOURCE_SP) <= DRIFT_SP + 1e-12
                assert se + sp - 1 >= ETA - 1e-12
        print(w.name, "true_delta", f"{delta(w.p_plus, w.p_minus):.12f}",
              "proxy_delta", f"{delta(w.q_plus, w.q_minus):.12f}")
    assert abs(delta(worlds[0].p_plus, worlds[0].p_minus) - 1 / 18) < 1e-12
    assert abs(delta(worlds[1].p_plus, worlds[1].p_minus) + 1 / 72) < 1e-12


def cp_interval(success: int, trials: int, alpha: float) -> tuple[float, float]:
    """Exact equal-tailed binomial interval, conditional on the observed denominator."""
    if trials == 0:
        return (0.0, 1.0)
    lo = 0.0 if success == 0 else float(beta.ppf(alpha / 2, success, trials - success + 1))
    hi = 1.0 if success == trials else float(beta.ppf(1 - alpha / 2, success + 1, trials - success))
    return (lo, hi)


def sample_primitives(rng: np.random.Generator, world: World, n_audit: int):
    gold_counts = [int(rng.binomial(N_SOURCE, p)) for p in SOURCE_P]
    source_p = [cp_interval(k, N_SOURCE, ALPHA / (7 if n_audit else 6)) for k in gold_counts]
    se_yes = sum(int(rng.binomial(k, SOURCE_SE)) for k in gold_counts)
    sp_yes = sum(int(rng.binomial(N_SOURCE - k, SOURCE_SP)) for k in gold_counts)
    k_primitives = 7 if n_audit else 6
    per_alpha = ALPHA / k_primitives
    source_se = cp_interval(se_yes, sum(gold_counts), per_alpha)
    source_sp = cp_interval(sp_yes, 2 * N_SOURCE - sum(gold_counts), per_alpha)
    q_plus = cp_interval(int(rng.binomial(N_TARGET, world.q_plus)), N_TARGET, per_alpha)
    q_minus = cp_interval(int(rng.binomial(N_TARGET, world.q_minus)), N_TARGET, per_alpha)
    audit = None
    if n_audit:
        audit = cp_interval(int(rng.binomial(n_audit, world.p_plus)), n_audit, per_alpha)
    return source_p, source_se, source_sp, (q_plus, q_minus), audit


def envelope(source_p, source_se, source_sp, q, audit=None):
    """Return conservative set projection [low, high] or None if infeasible.

    LP variables x=(u,v,p_plus,p_minus), with u=1/(Se+Sp-1), v=(1-Sp)u.
    All constraints are linear in x. For uncertain q, q_lo*u <= p+v <= q_hi*u.
    """
    s_l = max(0.0, source_se[0] - DRIFT_SE)
    s_u = min(1.0, source_se[1] + DRIFT_SE)
    c_l = max(0.0, source_sp[0] - DRIFT_SP)
    c_u = min(1.0, source_sp[1] + DRIFT_SP)
    p_bounds = [(max(0.0, ci[0] - DRIFT_P), min(1.0, ci[1] + DRIFT_P)) for ci in source_p]
    if audit is not None:
        p_bounds[0] = (max(p_bounds[0][0], audit[0]), min(p_bounds[0][1], audit[1]))
    if any(lo > hi for lo, hi in p_bounds):
        return None

    rows = [
        [s_l, -1, 0, 0],       # s_l*u-v <= 1
        [-s_u, 1, 0, 0],       # -s_u*u+v <= -1
        [1 - c_u, -1, 0, 0],   # (1-c_u)u <= v
        [-(1 - c_l), 1, 0, 0], # v <= (1-c_l)u
    ]
    rhs = [1, -1, 0, 0]
    for g, (q_lo, q_hi) in enumerate(q):
        row_lo = [-q_lo, 1, 0, 0]   # q_lo*u <= p+v
        row_hi = [q_hi, -1, 0, 0]   # p+v <= q_hi*u
        row_lo[2 + g] = 1
        row_hi[2 + g] = -1
        # row_lo represents -q_lo*u + v + p >= 0; flip for <=.
        rows.extend(([-x for x in row_lo], [-x for x in row_hi]))
        rhs.extend((0, 0))

    wp, wm = WEIGHTS
    c = np.array([0.0, 0.0, wp / (1 - T), -wm / (1 - T)])
    bounds = [(1.0, 1.0 / ETA), (None, None), p_bounds[0], p_bounds[1]]
    r_lo = linprog(c, A_ub=rows, b_ub=rhs, bounds=bounds, method="highs")
    if r_lo.status == 2:
        return None
    if not r_lo.success:
        raise RuntimeError(f"lower LP status={r_lo.status}: {r_lo.message}")
    r_hi = linprog(-c, A_ub=rows, b_ub=rhs, bounds=bounds, method="highs")
    if not r_hi.success:
        raise RuntimeError(f"upper LP status={r_hi.status}: {r_hi.message}")
    offset = -T * (wp - wm) / (1 - T)
    return (float(r_lo.fun + offset), float(-r_hi.fun + offset))


def exact_population_interval(world: World, n_audit: int):
    source_p = tuple((p, p) for p in SOURCE_P)
    audit = (world.p_plus, world.p_plus) if n_audit else None
    return envelope(source_p, (SOURCE_SE, SOURCE_SE), (SOURCE_SP, SOURCE_SP),
                    ((world.q_plus, world.q_plus), (world.q_minus, world.q_minus)), audit)


def run(n_reps: int, suffix: str) -> None:
    worlds = make_worlds()
    validate_worlds(worlds)
    pop_bounds = {w.name: {str(n): exact_population_interval(w, n) for n in AUDIT_SIZES}
                  for w in worlds}
    print("exact_population_intervals", pop_bounds, flush=True)
    if any(pop_bounds[w.name][str(n)] is None for w in worlds if w.assumed_valid for n in AUDIT_SIZES):
        raise AssertionError("Expected population model to be feasible")
    out = HERE / f"raw_replicates_{suffix}.csv"
    summary_out = HERE / f"summary_{suffix}.json"
    seeds = np.random.SeedSequence(SEED).spawn(len(worlds) * len(AUDIT_SIZES))
    summary = {
        "protocol": "ADEMP_LOCK.md", "seed": SEED, "n_repetitions_per_cell": n_reps,
        "n_source_per_group": N_SOURCE, "n_target_per_group": N_TARGET,
        "alpha": ALPHA, "population_intervals": pop_bounds, "cells": [],
    }
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=(
            "world", "n_audit", "rep", "lower", "upper", "true_delta", "true_covered",
            "population_set_covered", "certification", "correct_sign", "wrong_sign",
            "infeasible", "interval_width"))
        writer.writeheader()
        seed_i = 0
        for w in worlds:
            true_delta = delta(w.p_plus, w.p_minus)
            for n_audit in AUDIT_SIZES:
                pop_interval = pop_bounds[w.name][str(n_audit)]
                rng = np.random.default_rng(seeds[seed_i])
                seed_i += 1
                counts = {k: 0 for k in ("true_covered", "population_set_covered", "correct_sign",
                                          "wrong_sign", "infeasible", "indeterminate")}
                widths = []
                for rep in range(n_reps):
                    primitives = sample_primitives(rng, w, n_audit)
                    interval = envelope(*primitives)
                    infeasible = interval is None
                    if infeasible:
                        low = high = float("nan")
                        certification = "infeasible"
                    else:
                        low, high = interval
                        widths.append(high - low)
                        certification = "positive" if low > 1e-12 else "negative" if high < -1e-12 else "indeterminate"
                    true_covered = not infeasible and low - 1e-9 <= true_delta <= high + 1e-9
                    set_covered = (None if not w.assumed_valid or pop_interval is None else
                                   (not infeasible and low - 1e-9 <= pop_interval[0] and
                                    high + 1e-9 >= pop_interval[1]))
                    sign_true = "positive" if true_delta > 0 else "negative" if true_delta < 0 else "tie"
                    correct = certification == sign_true
                    wrong = certification in ("positive", "negative") and certification != sign_true
                    flags = dict(true_covered=true_covered, population_set_covered=set_covered,
                                 correct_sign=correct, wrong_sign=wrong, infeasible=infeasible,
                                 indeterminate=certification == "indeterminate")
                    for k, flag in flags.items():
                        if flag is not None:
                            counts[k] += int(flag)
                    writer.writerow(dict(world=w.name, n_audit=n_audit, rep=rep,
                                         lower=low, upper=high, true_delta=true_delta,
                                         true_covered=int(true_covered),
                                         population_set_covered=("" if set_covered is None else int(set_covered)),
                                         certification=certification, correct_sign=int(correct),
                                         wrong_sign=int(wrong), infeasible=int(infeasible),
                                         interval_width=(high - low if not infeasible else float("nan"))))
                metrics = {k: ({"proportion": v / n_reps,
                                "mcse": sqrt((v / n_reps) * (1 - v / n_reps) / n_reps),
                                "count": v} if not (k == "population_set_covered" and
                                                     (not w.assumed_valid or pop_interval is None))
                               else None) for k, v in counts.items()}
                metrics["mean_interval_width"] = float(np.mean(widths)) if widths else None
                metrics["median_interval_width"] = float(np.median(widths)) if widths else None
                cell = {"world": w.name, "assumptions_valid": w.assumed_valid,
                        "n_audit": n_audit, "true_delta": true_delta, "metrics": metrics}
                summary["cells"].append(cell)
                print(w.name, "audit", n_audit, "coverage", metrics["true_covered"]["proportion"],
                      "set_coverage", (metrics["population_set_covered"] or {}).get("proportion"),
                      "correct", metrics["correct_sign"]["proportion"],
                      "wrong", metrics["wrong_sign"]["proportion"], flush=True)
                summary_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("wrote", out, summary_out, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument("--suffix", required=True)
    args = parser.parse_args()
    if args.reps <= 0 or not args.suffix.replace("_", "").isalnum():
        parser.error("Positive --reps and alphanumeric/underscore --suffix required")
    run(args.reps, args.suffix)

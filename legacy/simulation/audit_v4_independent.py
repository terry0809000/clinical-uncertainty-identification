"""Independent row-level V4 audit from saved counts, without simulation imports.

Reconstructs seven CP intervals, both 4-variable LPs and all saved results.
Also verifies exact V3/V4 row correspondence and immutable V3 hashes.
"""

from __future__ import annotations

from collections import defaultdict
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from scipy.optimize import linprog
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
N = 10000
ALPHA = 0.05 / 7
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
METRICS = (
    "true_covered", "correct_model_set_covered", "correct_sign",
    "wrong_sign", "indeterminate", "infeasible",
)
WORLDS = {
    "B_negative": {
        "p": F(19, 60), "se": F(72, 100),
        "delta": F(-1, 72), "set": (F(-1, 72), F(-1, 72)),
    },
    "C_interior_tie": {
        "p": F(7, 20), "se": F(7, 10),
        "delta": F(0), "set": (F(-1, 152), F(1, 328)),
    },
}
EXPECTED_V3_HASHES = {
    "run_supplement_v3.py":
        "7B2B0496784018F1527F79AB8ECA8022DA6304972BCDF1E48A9703450628DFD6",
    "supplement_v3_raw_main_10000.csv":
        "9C04CF870C786A7C5798B148C3E8C8B260A595A4451D72CC87117A4A4D9660F3",
    "supplement_v3_summary_main_10000.json":
        "3649BC584F1E0E5AB33DC3B3E72DB6B1DCD78A01B07CA394DF681C67DEE6C8EB",
}


def close(a, b, tolerance=1e-12):
    return math.isclose(a, b, rel_tol=0, abs_tol=tolerance)


def cp(x: int, n: int):
    assert 0 <= x <= n
    if n == 0:
        return 0.0, 1.0
    lo = 0.0 if x == 0 else beta.ppf(ALPHA / 2, x, n - x + 1)
    hi = 1.0 if x == n else beta.ppf(1 - ALPHA / 2, x + 1, n - x)
    return float(lo), float(hi)


def reconstructed_envelope(bounds):
    """Fresh LP in variable order (p+, p-, inverse J, false-positive/J)."""
    pplus, pminus, se, sp, qplus, qminus, audit = bounds
    se_low = max(0.0, se[0] - .05)
    se_high = min(1.0, se[1] + .05)
    sp_low = max(0.0, sp[0] - .05)
    sp_high = min(1.0, sp[1] + .05)
    pplus_low = max(0.0, pplus[0] - .10, audit[0])
    pplus_high = min(1.0, pplus[1] + .10, audit[1])
    pminus_low = max(0.0, pminus[0] - .10)
    pminus_high = min(1.0, pminus[1] + .10)
    if pplus_low > pplus_high or pminus_low > pminus_high:
        return None
    # Se bounds: se_low*h-z <= 1; -se_high*h+z <= -1.
    # Sp bounds: (1-sp_high)h-z <=0; -(1-sp_low)h+z <=0.
    a = [
        [0, 0, se_low, -1],
        [0, 0, -se_high, 1],
        [0, 0, 1 - sp_high, -1],
        [0, 0, -(1 - sp_low), 1],
    ]
    b = [1, -1, 0, 0]
    for i, (lo, hi) in enumerate((qplus, qminus)):
        x = [0, 0, lo, -1]
        x[i] = -1
        a.append(x)  # q_low*h <= p_g+z
        x = [0, 0, -hi, 1]
        x[i] = 1
        a.append(x)  # p_g+z <= q_high*h
        b.extend((0, 0))
    limits = [
        (pplus_low, pplus_high), (pminus_low, pminus_high),
        (1.0, 2.0), (None, None),
    ]
    objective = [5 / 6, -5 / 12, 0, 0]
    low = linprog(objective, A_ub=a, b_ub=b, bounds=limits, method="highs")
    if low.status == 2:
        return None
    assert low.success, (low.status, low.message)
    high = linprog([-v for v in objective], A_ub=a, b_ub=b,
                   bounds=limits, method="highs")
    assert high.success, (high.status, high.message)
    return float(low.fun - 1 / 12), float(-high.fun - 1 / 12)


def expected_mechanism(world, lam, rse, rsp):
    p = WORLDS[world]["p"]
    se = WORLDS[world]["se"]
    lam = F(str(lam))
    rse = F(str(rse))
    rsp = F(str(rsp))
    q = F(31, 100)
    selected = (p + (lam - 1) * p * se) / (1 + (lam - 1) * q)
    recorded = rse * selected + (1 - rsp) * (1 - selected)
    return float(selected), float(recorded)


def main():
    for name, expected in EXPECTED_V3_HASHES.items():
        actual = hashlib.sha256((HERE / name).read_bytes()).hexdigest().upper()
        assert actual == expected, (name, actual, expected)

    v4_raw = HERE / "supplement_v4_raw_main_10000.csv"
    v3_raw = HERE / "supplement_v3_raw_main_10000.csv"
    v4_summary_file = HERE / "supplement_v4_summary_main_10000.json"
    v3_summary_file = HERE / "supplement_v3_summary_main_10000.json"
    v4_summary = json.loads(v4_summary_file.read_text(encoding="utf-8"))
    v3_summary = json.loads(v3_summary_file.read_text(encoding="utf-8"))
    assert v4_summary["protocol"] == "SUPPLEMENT_PROTOCOL_V4.md"
    assert v4_summary["seed"] == v3_summary["seed"] == 20260926
    assert v4_summary["repetitions_per_cell"] == v3_summary["repetitions_per_cell"] == N
    assert v4_summary["n_target_reference_audit"] == 800
    assert len(v4_summary["cells"]) == len(v3_summary["cells"]) == 8
    # The two summaries are otherwise identical by exact RNG replay.
    v4_without_protocol = {k: v for k, v in v4_summary.items() if k != "protocol"}
    v3_without_protocol = {k: v for k, v in v3_summary.items() if k != "protocol"}
    assert v4_without_protocol == v3_without_protocol

    counts = defaultdict(lambda: defaultdict(int))
    width_sums = defaultdict(float)
    width_ns = defaultdict(int)
    seen = defaultdict(set)
    cells = {
        (c["world"], c["lambda_proxy_positive"], c["reference_se"],
         c["reference_sp"]): c for c in v4_summary["cells"]
    }
    assert len(cells) == 8
    max_cp_diff = 0.0
    max_lp_diff = 0.0
    rows = 0
    with v4_raw.open(newline="", encoding="utf-8") as f4, \
            v3_raw.open(newline="", encoding="utf-8") as f3:
        r4 = csv.DictReader(f4)
        r3 = csv.DictReader(f3)
        assert tuple(r3.fieldnames) == BASE_FIELDS
        assert tuple(r4.fieldnames[:len(BASE_FIELDS)]) == BASE_FIELDS
        assert len(r4.fieldnames) == len(BASE_FIELDS) + 4 * len(PRIMITIVES)
        for row4, row3 in zip(r4, r3, strict=True):
            rows += 1
            for field in BASE_FIELDS:
                assert row4[field] == row3[field], (rows, field)
            key = (
                row4["world"], float(row4["lambda_proxy_positive"]),
                float(row4["ref_se"]), float(row4["ref_sp"]),
            )
            assert key in cells, key
            rep = int(row4["rep"])
            assert rep not in seen[key]
            seen[key].add(rep)
            assert key[0] in WORLDS
            exp_selected, exp_recorded = expected_mechanism(*key)
            assert close(float(row4["audit_selected_true_risk"]), exp_selected)
            assert close(float(row4["audit_recorded_positive_rate"]), exp_recorded)
            assert close(float(row4["true_delta"]), float(WORLDS[key[0]]["delta"]))
            obs = []
            for name in PRIMITIVES:
                x = int(row4[f"{name}_success"])
                n = int(row4[f"{name}_denominator"])
                lo, hi = cp(x, n)
                saved_lo = float(row4[f"{name}_lower"])
                saved_hi = float(row4[f"{name}_upper"])
                max_cp_diff = max(max_cp_diff, abs(saved_lo - lo), abs(saved_hi - hi))
                assert close(saved_lo, lo, 2e-14), (rows, name, "lo")
                assert close(saved_hi, hi, 2e-14), (rows, name, "hi")
                obs.append((lo, hi))
            assert [int(row4[f"{name}_denominator"]) for name in PRIMITIVES] == [
                1000, 1000,
                int(row4["source_gold_plus_success"]) +
                int(row4["source_gold_minus_success"]),
                2000 - int(row4["source_gold_plus_success"]) -
                int(row4["source_gold_minus_success"]),
                1000, 1000, 800,
            ]
            recomputed = reconstructed_envelope(obs)
            infeasible = recomputed is None
            assert infeasible == bool(int(row4["infeasible"])), (rows, "feasibility")
            if infeasible:
                assert all(math.isnan(float(row4[x])) for x in ("lower", "upper", "width"))
                certification = "infeasible"
            else:
                lo, hi = recomputed
                saved_lo, saved_hi = float(row4["lower"]), float(row4["upper"])
                max_lp_diff = max(max_lp_diff, abs(lo - saved_lo), abs(hi - saved_hi))
                assert close(lo, saved_lo, 1e-8), (rows, "LP lower", lo, saved_lo)
                assert close(hi, saved_hi, 1e-8), (rows, "LP upper", hi, saved_hi)
                assert close(float(row4["width"]), hi - lo, 1e-8)
                width_sums[key] += float(row4["width"])
                width_ns[key] += 1
                # Certification is checked against the saved output and LP endpoints.
                certification = (
                    "positive" if saved_lo > 1e-12 else
                    "negative" if saved_hi < -1e-12 else "indeterminate"
                )
            assert row4["certification"] == certification
            truth = float(WORLDS[key[0]]["delta"])
            population_set = WORLDS[key[0]]["set"]
            saved_lo = float(row4["lower"])
            saved_hi = float(row4["upper"])
            true_sign = (
                "tie" if truth == 0 else "positive" if truth > 0 else "negative"
            )
            flags = {
                "true_covered": (
                    not infeasible and saved_lo - 1e-9 <= truth <= saved_hi + 1e-9
                ),
                "correct_model_set_covered": (
                    not infeasible
                    and saved_lo - 1e-9 <= float(population_set[0])
                    and saved_hi + 1e-9 >= float(population_set[1])
                ),
                "correct_sign": certification == true_sign,
                "wrong_sign": (
                    certification in ("positive", "negative")
                    and certification != true_sign
                ),
                "indeterminate": certification == "indeterminate",
                "infeasible": infeasible,
            }
            for metric, flag in flags.items():
                assert int(row4[metric]) == int(flag), (rows, metric)
                counts[key][metric] += int(flag)
    assert rows == 8 * N
    assert len(seen) == 8 and all(reps == set(range(N)) for reps in seen.values())
    for key, cell in cells.items():
        assert sum(counts[key][m] for m in
                   ("correct_sign", "wrong_sign", "indeterminate", "infeasible")) == N
        for metric in METRICS:
            x = counts[key][metric]
            item = cell["metrics"][metric]
            p = x / N
            assert item["count"] == x, (key, metric)
            assert close(item["proportion"], p, 1e-15)
            assert close(item["mcse"], math.sqrt(p * (1 - p) / N), 1e-15)
            expected_upper = (
                1.0 if x == N else float(beta.ppf(.95, x + 1, N - x))
            )
            assert close(item["one_sided_95_upper"], expected_upper, 1e-14)
        assert close(cell["metrics"]["mean_width_when_feasible"],
                     width_sums[key] / width_ns[key], 1e-12)
        exp_selected, exp_recorded = expected_mechanism(*key)
        assert close(cell["selected_true_risk"], exp_selected)
        assert close(cell["recorded_audit_positive_rate"], exp_recorded)
        assert close(cell["true_delta"], float(WORLDS[key[0]]["delta"]))
        for a, b in zip(cell["correct_model_population_set"],
                        WORLDS[key[0]]["set"]):
            assert close(a, float(b), 1e-12)
        assert cell["assumptions_valid"] == (key[1:] == (1.0, 1.0, 1.0))
    assert len(v4_summary["interactions"]) == 4
    for item in v4_summary["interactions"]:
        world, metric = item["world"], item["outcome"]
        assert metric in ("true_covered", "wrong_sign")
        rates = {}
        var = {}
        for key in cells:
            if key[0] != world:
                continue
            factor = (key[1] == 4.0, key[2] < 1.0)
            p = counts[key][metric] / N
            rates[factor] = p
            var[factor] = p * (1 - p) / N
        assert len(rates) == 4
        interaction = (
            rates[True, True] - rates[True, False]
            - rates[False, True] + rates[False, False]
        )
        se = math.sqrt(sum(var.values()))
        assert close(item["probability_scale_interaction"], interaction, 1e-14)
        assert close(item["mcse"], se, 1e-14)
        assert close(item["descriptive_95_monte_carlo_interval"][0],
                     interaction - 1.96 * se, 1e-14)
        assert close(item["descriptive_95_monte_carlo_interval"][1],
                     interaction + 1.96 * se, 1e-14)

    print("RAW_ROWS", rows)
    print("MAX_CP_ABSOLUTE_DIFFERENCE", format(max_cp_diff, ".17g"))
    print("MAX_LP_ABSOLUTE_DIFFERENCE", format(max_lp_diff, ".17g"))
    for key in sorted(cells):
        print("CELL", key, "covered", counts[key]["true_covered"],
              "wrong", counts[key]["wrong_sign"],
              "infeasible", counts[key]["infeasible"])
    for path in (
        HERE / "SUPPLEMENT_PROTOCOL_V4.md",
        HERE / "run_supplement_v4.py",
        HERE / "audit_v4_independent.py",
        v4_raw, v4_summary_file,
    ):
        print("SHA256", path.name, hashlib.sha256(path.read_bytes()).hexdigest().upper())
    print("ALL_V4_ROW_LEVEL_CHECKS_PASSED")


if __name__ == "__main__":
    main()

"""Independent, standard-library audit of the V3 raw Monte Carlo record.

Does not import the simulation or projection implementation. Raises if any
record-level flag or aggregated metric conflicts with the saved summary.
"""

from collections import defaultdict
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
RAW = HERE / "supplement_v3_raw_main_10000.csv"
SUMMARY = HERE / "supplement_v3_summary_main_10000.json"
N = 10000
WORLDS = {
    "B_negative": {
        "p_plus": F(19, 60), "se": F(72, 100), "sp": F(88, 100),
        "delta": F(-1, 72), "set": (F(-1, 72), F(-1, 72)),
    },
    "C_interior_tie": {
        "p_plus": F(7, 20), "se": F(7, 10), "sp": F(9, 10),
        "delta": F(0), "set": (F(-1, 152), F(1, 328)),
    },
}
METRICS = (
    "true_covered", "correct_model_set_covered", "correct_sign",
    "wrong_sign", "indeterminate", "infeasible",
)


def cell_key(world, lam, ref_se, ref_sp):
    return world, float(lam), float(ref_se), float(ref_sp)


def mechanism(w, lam, ref_se, ref_sp):
    p, se, sp = w["p_plus"], w["se"], w["sp"]
    q = p * se + (1 - p) * (1 - sp)
    assert q == F(31, 100)
    # Independent derivation via selection-weighted joint true-positive mass.
    selected_p = (p + (lam - 1) * p * se) / (1 + (lam - 1) * q)
    recorded_p = ref_se * selected_p + (1 - ref_sp) * (1 - selected_p)
    return float(selected_p), float(recorded_p)


summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
assert summary["seed"] == 20260926
assert summary["repetitions_per_cell"] == N
assert summary["n_target_reference_audit"] == 800
assert len(summary["cells"]) == 8
cells = {}
for c in summary["cells"]:
    key = cell_key(c["world"], c["lambda_proxy_positive"],
                   c["reference_se"], c["reference_sp"])
    assert key not in cells
    cells[key] = c
    assert c["assumptions_valid"] == (
        key[1:] == (1.0, 1.0, 1.0)
    )
    w = WORLDS[key[0]]
    expected_selected, expected_recorded = mechanism(
        w, F(str(key[1])), F(str(key[2])), F(str(key[3]))
    )
    assert math.isclose(c["selected_true_risk"], expected_selected, abs_tol=1e-14)
    assert math.isclose(c["recorded_audit_positive_rate"], expected_recorded, abs_tol=1e-14)
    assert math.isclose(c["true_delta"], float(w["delta"]), abs_tol=1e-14)
    for observed, expected in zip(c["correct_model_population_set"], w["set"]):
        assert math.isclose(observed, float(expected), abs_tol=1e-12)

seen = defaultdict(set)
counts = defaultdict(lambda: defaultdict(int))
width_sum = defaultdict(float)
width_count = defaultdict(int)
line_count = 0
with RAW.open("r", encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    assert reader.fieldnames is not None and len(reader.fieldnames) == 18
    for row in reader:
        line_count += 1
        key = cell_key(row["world"], row["lambda_proxy_positive"],
                       row["ref_se"], row["ref_sp"])
        assert key in cells, key
        rep = int(row["rep"])
        assert rep not in seen[key], (key, rep)
        seen[key].add(rep)
        c = cells[key]
        w = WORLDS[key[0]]
        expected_selected, expected_recorded = mechanism(
            w, F(str(key[1])), F(str(key[2])), F(str(key[3]))
        )
        assert math.isclose(float(row["audit_selected_true_risk"]),
                            expected_selected, abs_tol=1e-14)
        assert math.isclose(float(row["audit_recorded_positive_rate"]),
                            expected_recorded, abs_tol=1e-14)
        truth = float(row["true_delta"])
        assert math.isclose(truth, float(w["delta"]), abs_tol=1e-14)
        low, high, width = (float(row[x]) for x in ("lower", "upper", "width"))
        infeasible = row["certification"] == "infeasible"
        if infeasible:
            assert math.isnan(low) and math.isnan(high) and math.isnan(width)
            certification = "infeasible"
        else:
            assert math.isfinite(low) and math.isfinite(high) and low <= high
            assert math.isclose(width, high - low, abs_tol=1e-13)
            width_sum[key] += width
            width_count[key] += 1
            certification = (
                "positive" if low > 1e-12 else
                "negative" if high < -1e-12 else "indeterminate"
            )
        assert row["certification"] == certification
        true_sign = "tie" if w["delta"] == 0 else (
            "positive" if w["delta"] > 0 else "negative"
        )
        expected_flags = {
            "true_covered": (not infeasible and low - 1e-9 <= truth <= high + 1e-9),
            "correct_model_set_covered": (
                not infeasible and low - 1e-9 <= float(w["set"][0])
                and high + 1e-9 >= float(w["set"][1])
            ),
            "correct_sign": certification == true_sign,
            "wrong_sign": (certification in ("positive", "negative")
                           and certification != true_sign),
            "indeterminate": certification == "indeterminate",
            "infeasible": infeasible,
        }
        for metric, flag in expected_flags.items():
            assert int(row[metric]) == int(flag), (key, rep, metric)
            counts[key][metric] += int(flag)
assert line_count == 8 * N
assert len(seen) == 8 and all(reps == set(range(N)) for reps in seen.values())

for key, c in cells.items():
    assert sum(counts[key][x] for x in
               ("correct_sign", "wrong_sign", "indeterminate", "infeasible")) == N
    for metric in METRICS:
        x = counts[key][metric]
        item = c["metrics"][metric]
        assert item["count"] == x, (key, metric, x, item["count"])
        assert math.isclose(item["proportion"], x / N, abs_tol=1e-15)
        assert math.isclose(item["mcse"], math.sqrt((x / N) * (1 - x / N) / N),
                            abs_tol=1e-15)
    assert math.isclose(c["metrics"]["mean_width_when_feasible"],
                        width_sum[key] / width_count[key], abs_tol=1e-12)

interaction_output = {}
for item in summary["interactions"]:
    world, metric = item["world"], item["outcome"]
    assert metric in ("true_covered", "wrong_sign")
    rates, variances = {}, {}
    for key in cells:
        if key[0] != world:
            continue
        ix = (key[1] == 4.0, key[2] == 0.8)
        p = counts[key][metric] / N
        rates[ix] = p
        variances[ix] = p * (1 - p) / N
    assert len(rates) == 4
    interaction = rates[True, True] - rates[True, False] - rates[False, True] + rates[False, False]
    mcse = math.sqrt(sum(variances.values()))
    assert math.isclose(item["probability_scale_interaction"], interaction, abs_tol=1e-14)
    assert math.isclose(item["mcse"], mcse, abs_tol=1e-14)
    for obs, exp in zip(item["descriptive_95_monte_carlo_interval"],
                        (interaction - 1.96 * mcse, interaction + 1.96 * mcse)):
        assert math.isclose(obs, exp, abs_tol=1e-14)
    interaction_output[f"{world}/{metric}"] = [interaction, mcse]

print("RAW_ROWS", line_count)
for key in sorted(cells):
    c = cells[key]
    print("CELL", key,
          "selected", format(c["selected_true_risk"], ".10f"),
          "recorded", format(c["recorded_audit_positive_rate"], ".10f"),
          "covered", counts[key]["true_covered"],
          "wrong", counts[key]["wrong_sign"],
          "infeasible", counts[key]["infeasible"],
          "whole_set", counts[key]["correct_model_set_covered"])
for key, pair in sorted(interaction_output.items()):
    print("INTERACTION", key, pair)
for p in (RAW, SUMMARY):
    print("SHA256", p.name, hashlib.sha256(p.read_bytes()).hexdigest())
print("ALL_V3_RAW_AUDIT_CHECKS_PASSED")

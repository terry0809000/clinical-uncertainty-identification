"""Exact-rational 56-cell sensitivity surface for Box 1.

Design locked in SENSITIVITY_SURFACE_PROTOCOL_V5.md before this code/results.
This script only manipulates stipulated synthetic population probabilities.
"""

from __future__ import annotations

import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
Q = (F(31, 100), F(2, 5))
P0 = (F(2, 5), F(11, 20))
SE0, SP0, ETA = F(67, 100), F(93, 100), F(1, 2)
DP_GRID = tuple(map(F, ("0", ".01", ".02", ".03", ".05", ".075", ".10", ".15")))
DC_GRID = tuple(map(F, ("0", ".005", ".01", ".02", ".03", ".04", ".05")))


def constraints(dp: F, dc: F) -> list[tuple[F, F, F, str]]:
    # All rows mean A*a + B*b <= C, a=Se+Sp-1, b=1-Sp.
    rows = [
        (F(-1), F(0), -ETA, "Youden >= 1/2"),
        (F(0), F(-1), -(1 - SP0 - dc), "Sp upper"),
        (F(0), F(1), 1 - SP0 + dc, "Sp lower"),
        (F(-1), F(-1), -(SE0 - dc), "Se lower"),
        (F(1), F(1), SE0 + dc, "Se upper"),
        (F(0), F(-1), F(0), "b >= 0"),
        (F(0), F(1), F(1), "b <= 1"),
        (F(-1), F(-1), F(0), "Se >= 0"),
        (F(1), F(1), F(1), "Se <= 1"),
    ]
    for q, p0, group in zip(Q, P0, ("+", "-")):
        lo, hi = max(F(0), p0 - dp), min(F(1), p0 + dp)
        rows.extend((
            (lo, F(1), q, f"p{group} >= {lo}"),
            (-hi, F(-1), -q, f"p{group} <= {hi}"),
        ))
    return rows


def pair_intersection(one, two):
    A, B, C, _ = one
    D, E, H, _ = two
    det = A * E - B * D
    if not det:
        return None
    return ((C * E - B * H) / det, (A * H - C * D) / det)


def feasible(point: tuple[F, F], rows) -> bool:
    a, b = point
    return a >= ETA and all(A * a + B * b <= C for A, B, C, _ in rows)


def world(a: F, b: F) -> dict:
    pplus, pminus = (Q[0] - b) / a, (Q[1] - b) / a
    se, sp = a + b, 1 - b
    delta = (10 * pplus - 5 * pminus - 1) / 12
    for p, q in ((pplus, Q[0]), (pminus, Q[1])):
        assert 0 <= p <= 1 and 0 <= se <= 1 and 0 <= sp <= 1
        joint = (p * se, p * (1 - se), (1 - p) * (1 - sp), (1 - p) * sp)
        assert all(x >= 0 for x in joint) and sum(joint) == 1
        assert joint[0] + joint[2] == q
    return {"a": a, "b": b, "pplus": pplus, "pminus": pminus,
            "se": se, "sp": sp, "delta": delta}


def serial(x: F) -> dict:
    return {"fraction": f"{x.numerator}/{x.denominator}", "decimal": float(x)}


def serial_world(d: dict) -> dict:
    return {k: serial(v) for k, v in d.items()}


def main() -> None:
    records = []
    for dp in DP_GRID:
        for dc in DC_GRID:
            rows = constraints(dp, dc)
            vertices = sorted({point for pair in combinations(rows, 2)
                               if (point := pair_intersection(*pair)) is not None
                               and feasible(point, rows)})
            assert vertices, (dp, dc)
            worlds = [world(*point) for point in vertices]
            low = min(worlds, key=lambda x: x["delta"])
            high = max(worlds, key=lambda x: x["delta"])
            lower, upper = low["delta"], high["delta"]
            status = ("strict_positive" if lower > 0 else
                      "strict_negative" if upper < 0 else "zero_compatible")
            assert lower <= F(1, 48) <= upper
            records.append({
                "dp": serial(dp), "dc": serial(dc), "vertex_count": len(vertices),
                "lower": serial(lower), "upper": serial(upper),
                "width": serial(upper - lower), "sign_status": status,
                "lower_world": serial_world(low), "upper_world": serial_world(high),
            })

    # Enlarging either sensitivity class must expand, never shrink, the set.
    lookup = {(F(r["dp"]["fraction"]), F(r["dc"]["fraction"])): r for r in records}
    for dp in DP_GRID:
        for dc in DC_GRID:
            cur = lookup[dp, dc]
            for next_dp in DP_GRID:
                if next_dp > dp:
                    enlarged = lookup[next_dp, dc]
                    assert F(enlarged["lower"]["fraction"]) <= F(cur["lower"]["fraction"])
                    assert F(enlarged["upper"]["fraction"]) >= F(cur["upper"]["fraction"])
            for next_dc in DC_GRID:
                if next_dc > dc:
                    enlarged = lookup[dp, next_dc]
                    assert F(enlarged["lower"]["fraction"]) <= F(cur["lower"]["fraction"])
                    assert F(enlarged["upper"]["fraction"]) >= F(cur["upper"]["fraction"])

    box = lookup[F(1, 10), F(1, 20)]
    assert (F(box["lower"]["fraction"]), F(box["upper"]["fraction"])) == (F(-1, 72), F(1, 18))
    origin = lookup[F(0), F(0)]
    assert F(origin["lower"]["fraction"]) == F(origin["upper"]["fraction"]) == F(1, 48)

    first_zero_by_dc = {}
    for dc in DC_GRID:
        first = next((dp for dp in DP_GRID if lookup[dp, dc]["sign_status"] == "zero_compatible"), None)
        first_zero_by_dc[str(dc)] = None if first is None else serial(first)
    result = {
        "status": "exact synthetic population sensitivity surface; independent check pending",
        "protocol": "SENSITIVITY_SURFACE_PROTOCOL_V5.md",
        "cells": len(records),
        "dp_grid": [serial(x) for x in DP_GRID],
        "dc_grid": [serial(x) for x in DC_GRID],
        "first_zero_compatible_dp_by_dc_on_grid": first_zero_by_dc,
        "rows": records,
    }
    dest = HERE / "sensitivity_surface_v5_exact.json"
    dest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(dest)
    print("cells", len(records), "sign_counts",
          {s: sum(x["sign_status"] == s for x in records)
           for s in ("strict_positive", "strict_negative", "zero_compatible")})
    print("first_zero_compatible_dp_by_dc_on_grid", first_zero_by_dc)
    print("box", box["lower"], box["upper"])


if __name__ == "__main__":
    main()

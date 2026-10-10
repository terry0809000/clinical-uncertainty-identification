"""Numerical cross-check of V5 exact vertices using a different LP parameterization.

This file does not import run_sensitivity_surface_v5.py or reuse its constraints.
Its result is a computational cross-check, not external scientific review.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


HERE = Path(__file__).resolve().parent
source = HERE / "sensitivity_surface_v5_exact.json"
saved = json.loads(source.read_text(encoding="utf-8"))


def solve(dp: float, dc: float, maximize: bool) -> float:
    # Charnes-Cooper variables u=1/a, v=b/a; target risks are q*u-v.
    # delta=(11*u-50*v-10)/120, a linear objective in (u,v).
    qplus, qminus = .31, .40
    p0plus, p0minus = .40, .55
    se0, sp0 = .67, .93
    se_lo, se_hi = max(0, se0 - dc), min(1, se0 + dc)
    sp_lo, sp_hi = max(0, sp0 - dc), min(1, sp0 + dc)
    rows = []
    rhs = []

    def add(coeff_u, coeff_v, bound):
        rows.append([coeff_u, coeff_v])
        rhs.append(bound)

    # v/u = b = 1-Sp, and Se = (1+v)/u.
    add(1 - sp_hi, -1, 0)  # v >= (1-sp_hi)*u
    add(-(1 - sp_lo), 1, 0)  # v <= (1-sp_lo)*u
    add(-se_hi, 1, -1)  # 1+v <= se_hi*u
    add(se_lo, -1, 1)  # 1+v >= se_lo*u
    add(-1, 1, 0)  # b <= 1

    for q, p0 in ((qplus, p0plus), (qminus, p0minus)):
        lo, hi = max(0, p0 - dp), min(1, p0 + dp)
        add(-q, 1, -lo)  # q*u-v >= lo
        add(q, -1, hi)  # q*u-v <= hi

    c = np.array([11 / 120, -50 / 120], dtype=float)
    if maximize:
        c = -c
    res = linprog(c, A_ub=rows, b_ub=rhs, bounds=[(1, 2), (0, None)], method="highs")
    if not res.success:
        raise RuntimeError((dp, dc, maximize, res.message))
    return float((11 * res.x[0] - 50 * res.x[1] - 10) / 120)


errors = []
for row in saved["rows"]:
    dp, dc = row["dp"]["decimal"], row["dc"]["decimal"]
    low, high = solve(dp, dc, False), solve(dp, dc, True)
    exact_low, exact_high = row["lower"]["decimal"], row["upper"]["decimal"]
    errors.append(max(abs(low - exact_low), abs(high - exact_high)))
max_error = max(errors)
assert len(errors) == 56 and max_error < 1e-10, (len(errors), max_error)

result = {
    "pass": True,
    "scope": "56-cell independent numerical Charnes-Cooper LP cross-check; not external scientific review",
    "exact_json_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "cells": len(errors),
    "maximum_absolute_endpoint_difference": max_error,
    "all_endpoint_differences_below_1e-10": True,
}
out = HERE / "sensitivity_surface_v5_lp_audit.json"
out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))

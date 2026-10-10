"""Exact rational audit of Box 1's deterministic identification example.

This enumerates vertices of a two-dimensional rational polytope and its
intersections with an exact one-stratum audit line. It does not analyse data.
"""

from fractions import Fraction as F
from itertools import combinations


Q_PLUS, Q_MINUS = F(31, 100), F(2, 5)
P0_PLUS, P0_MINUS = F(2, 5), F(11, 20)
DP, DS, DC = F(1, 10), F(1, 20), F(1, 20)
SE0, SP0, ETA = F(67, 100), F(93, 100), F(1, 2)

# Each inequality is A*a + B*b <= C.
constraints = [
    (-1, 0, -ETA, "a >= eta"),
    (0, -1, -(1 - SP0 - DC), "b >= 1-Sp_upper"),
    (0, 1, 1 - SP0 + DC, "b <= 1-Sp_lower"),
    (-1, -1, -(SE0 - DS), "Se >= Se_lower"),
    (1, 1, SE0 + DS, "Se <= Se_upper"),
]
for q, p0, label in ((Q_PLUS, P0_PLUS, "+"),
                     (Q_MINUS, P0_MINUS, "-")):
    lo, hi = p0 - DP, p0 + DP
    constraints.extend([
        (lo, 1, q, f"p{label} >= {lo}"),
        (-hi, -1, -q, f"p{label} <= {hi}"),
    ])


def feasible(a, b):
    return all(F(A) * a + F(B) * b <= F(C)
               for A, B, C, _ in constraints)


def vertex(eq1, eq2):
    A, B, C, _ = eq1
    D, E, H, _ = eq2
    det = F(A) * E - F(B) * D
    if det == 0:
        return None
    return ((F(C) * E - F(B) * H) / det,
            (F(A) * H - F(C) * D) / det)


def values(a, b):
    pp, pm = (Q_PLUS - b) / a, (Q_MINUS - b) / a
    se, sp = a + b, 1 - b
    delta = (10 * pp - 5 * pm - 1) / 12
    return pp, pm, se, sp, delta


vertices = sorted({point for pair in combinations(constraints, 2)
                   if (point := vertex(*pair)) is not None and feasible(*point)})
assert vertices
print("NO_AUDIT_VERTICES")
for a, b in vertices:
    print("a", a, "b", b, "(p+,p-,Se,Sp,Delta)", values(a, b))

no_audit = (min(values(*v)[-1] for v in vertices),
            max(values(*v)[-1] for v in vertices))
print("NO_AUDIT_SET", no_audit)
assert no_audit == (F(-1, 72), F(1, 18))


def exact_audit(rate):
    # b = q+ - a*rate. A linear-fractional objective on this line
    # is monotone in a here, so endpoints occur at constraint crossings.
    candidate = set()
    for A, B, C, _ in constraints:
        denom = F(A) - F(B) * rate
        if denom:
            a = (F(C) - F(B) * Q_PLUS) / denom
            b = Q_PLUS - a * rate
            if a > 0 and feasible(a, b):
                candidate.add((a, b))
    candidate = sorted(candidate)
    assert candidate
    for a, b in candidate:
        assert (Q_PLUS - b) / a == rate
    bounds = (min(values(*x)[-1] for x in candidate),
              max(values(*x)[-1] for x in candidate))
    return candidate, bounds


worlds = {
    "A_positive": F(29, 60),
    "B_negative": F(19, 60),
    "C_interior_tie": F(7, 20),
}
expected = {
    "A_positive": (F(1, 18), F(1, 18)),
    "B_negative": (F(-1, 72), F(-1, 72)),
    "C_interior_tie": (F(-1, 152), F(1, 328)),
}
for label, rate in worlds.items():
    candidate, bounds = exact_audit(rate)
    print(label, "p+", rate, "SET", bounds)
    for a, b in candidate:
        print("  endpoint a", a, "b", b,
              "(p+,p-,Se,Sp,Delta)", values(a, b))
    assert bounds == expected[label], (label, bounds)

print("ALL_EXACT_CHECKS_PASSED")

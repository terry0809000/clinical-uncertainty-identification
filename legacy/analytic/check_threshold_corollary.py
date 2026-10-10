"""Exact check of the post hoc fixed-rule threshold corollary.

This reuses only the already established Box 1 sharp endpoints. It does not
claim an independent proof of the original feasible-set optimization.
"""

from fractions import Fraction as F


def net_benefit(p_plus: F, p_minus: F, threshold: F) -> F:
    assert F(0) < threshold < F(1)
    return (2 * p_plus - p_minus - threshold) / (3 * (1 - threshold))


def interval(threshold: F) -> tuple[F, F]:
    assert F(0) < threshold < F(1)
    return ((F(1, 6) - threshold) / (3 * (1 - threshold)),
            (F(1, 3) - threshold) / (3 * (1 - threshold)))


def main() -> None:
    worlds = {
        "A": (F(29, 60), F(19, 30)),
        "B": (F(19, 60), F(7, 15)),
        "C": (F(7, 20), F(1, 2)),
    }
    score = {name: 2 * pp - pm for name, (pp, pm) in worlds.items()}
    assert score == {"A": F(1, 3), "B": F(1, 6), "C": F(1, 5)}
    assert interval(F(1, 5)) == (F(-1, 72), F(1, 18))
    baseline_low, baseline_high = interval(F(1, 5))
    worst_regret_rule_1 = max(F(0), -baseline_low)
    worst_regret_rule_2 = max(F(0), baseline_high)
    assert (worst_regret_rule_1, worst_regret_rule_2) == (F(1, 72), F(1, 18))
    assert worst_regret_rule_1 < worst_regret_rule_2
    assert (max(F(0), -interval(F(1, 4))[0]),
            max(F(0), interval(F(1, 4))[1])) == (F(1, 27), F(1, 27))
    for t in (F(1, 10), F(1, 6), F(1, 5), F(1, 4), F(1, 3), F(1, 2)):
        lower, upper = interval(t)
        assert lower == net_benefit(*worlds["B"], t)
        assert upper == net_benefit(*worlds["A"], t)
        assert lower <= net_benefit(*worlds["C"], t) <= upper
        expected = (
            "positive" if t < F(1, 6) else
            "nonnegative_with_tie" if t == F(1, 6) else
            "both_signs" if t < F(1, 3) else
            "nonpositive_with_tie" if t == F(1, 3) else
            "negative"
        )
        actual = (
            "positive" if lower > 0 else
            "negative" if upper < 0 else
            "nonnegative_with_tie" if lower == 0 else
            "nonpositive_with_tie" if upper == 0 else
            "both_signs"
        )
        assert actual == expected, (t, actual, expected)
        print(f"t={t}: [{lower},{upper}], {actual}")
    print("PASS: affine mapping and sign boundaries checked with exact fractions")


if __name__ == "__main__":
    main()

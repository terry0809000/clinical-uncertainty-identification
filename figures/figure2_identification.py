"""Figure 2: exact analytic identified sets for observationally equivalent worlds.

The worlds are synthetic mathematical constructions. No clinical outcome data are
used. Endpoints are computed by the independently audited LP implementation.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from sys import path as syspath
syspath.insert(0, str(Path(__file__).resolve().parents[1] / "legacy" / "simulation"))
from run_study import delta, exact_population_interval, make_worlds
from run_supplement_v2 import new_worlds

OUT = Path(__file__).resolve().parent
A, B, _ = make_worlds()
C, _ = new_worlds()
worlds = (A, B, C)
for w in worlds:
    assert abs(w.q_plus - 0.31) < 1e-12 and abs(w.q_minus - 0.40) < 1e-12

COL = {"navy": "#13253D", "blue": "#176BB0", "coral": "#B83B4D",
       "violet": "#6947A8", "teal": "#08756F", "muted": "#5C6B7B",
       "line": "#D9E2EB", "paper": "#FFFFFF"}
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.spines.right": False,
                     "axes.spines.top": False, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})

fig = plt.figure(figsize=(6.69, 3.6), facecolor=COL["paper"])
axq = fig.add_axes((0.10, 0.30, 0.29, 0.52))
axd = fig.add_axes((0.51, 0.30, 0.43, 0.52))

# A: identical observable target proxy evidence for all three worlds.
axq.bar([0, 1], [0.31, 0.40], width=0.5,
        color=[COL["blue"], COL["teal"]], edgecolor=COL["navy"], linewidth=0.7)
for x, y in ((0, 0.31), (1, 0.40)):
    axq.text(x, y + 0.012, f"{y:.2f}", ha="center", va="bottom", fontsize=8,
             fontweight="bold", color=COL["navy"])
axq.set_xticks([0, 1], [r"$G_+$", r"$G_-$"])
axq.set_ylim(0, 0.49)
axq.set_ylabel("Recorded outcome rate")
axq.set_title("A   Same observed target labels", loc="left", color=COL["navy"],
              fontsize=9.1, fontweight="bold", pad=7)
axq.spines["left"].set_color(COL["line"])
axq.spines["bottom"].set_color(COL["line"])
axq.tick_params(axis="both", length=0, labelcolor=COL["muted"])
axq.set_yticks([0, 0.2, 0.4])
axq.grid(axis="y", color=COL["line"], linewidth=0.55, zorder=0)
axq.set_axisbelow(True)
axq.text(0.5, -0.28, "Worlds A, B and C also share\nsource audited probabilities.",
         transform=axq.transAxes, ha="center", va="top", color=COL["muted"], fontsize=7.1)

# B: an exact sharp interval before audit and distinct post-audit sets.
un = exact_population_interval(A, 0)
assert all(abs(x-y) < 1e-12 for w in (B, C) for x, y in zip(exact_population_interval(w, 0), un))
rows = [
    ("No target audit", un, COL["navy"], None),
    ("Audit A", exact_population_interval(A, 50), COL["teal"], delta(A.p_plus, A.p_minus)),
    ("Audit C", exact_population_interval(C, 50), COL["violet"], delta(C.p_plus, C.p_minus)),
    ("Audit B", exact_population_interval(B, 50), COL["coral"], delta(B.p_plus, B.p_minus)),
]
for yy, (label, interval, color, point) in zip([3, 2, 1, 0], rows):
    lo, hi = interval
    axd.plot([lo, hi], [yy, yy], color=color, linewidth=5.0,
             solid_capstyle="round", zorder=3)
    if point is not None:
        axd.scatter([point], [yy], marker="o" if label != "Audit C" else "D",
                    s=38, color=color, edgecolor="white", linewidth=0.7, zorder=4)
axd.axvline(0, color=COL["muted"], linewidth=0.9, linestyle=(0, (3, 2)), zorder=1)
axd.set_yticks([3, 2, 1, 0], [r[0] for r in rows])
axd.set_xlim(-0.024, 0.067)
axd.set_xticks([-0.02, 0.00, 0.02, 0.04, 0.06])
axd.set_xlabel(r"Target net-benefit difference  $\Delta NB$")
axd.set_ylim(-0.5, 3.5)
axd.set_title("B   Identification after reference audit", loc="left",
              color=COL["navy"], fontsize=9.1, fontweight="bold", pad=7)
axd.spines["left"].set_visible(False)
axd.spines["bottom"].set_color(COL["line"])
axd.tick_params(axis="y", length=0, labelcolor=COL["navy"], labelsize=7.4)
axd.tick_params(axis="x", colors=COL["muted"], labelsize=7)
axd.grid(axis="x", color=COL["line"], linewidth=0.55, zorder=0)
axd.set_axisbelow(True)
axd.text(0.50, -0.28, "Exact population sets; audit rows use\nworld-specific exact $p_{1+}$ (not sample estimates).",
         transform=axd.transAxes, ha="center", va="top", color=COL["muted"], fontsize=7.1)

# Article-like figure title is kept in the manuscript legend, not the artwork.
for ext in ("pdf", "svg", "eps", "png"):
    fig.savefig(OUT / f"figure_2.{ext}", dpi=600 if ext == "png" else None)
plt.close(fig)
print("Analytic unaudited set", un)
for w in worlds:
    print(w.name, "delta", delta(w.p_plus, w.p_minus),
          "audit set", exact_population_interval(w, 50))

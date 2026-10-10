"""Post hoc exact-algebra figure for the fixed-rule threshold corollary.

This visualizes stipulated synthetic population bounds and regret, not patient
results or a newly selected clinical threshold. Rules and feasible laws remain
fixed while the utility threshold changes.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parent
T = np.linspace(0.05, 0.50, 700)
LOW = (1 / 6 - T) / (3 * (1 - T))
HIGH = (1 / 3 - T) / (3 * (1 - T))
REGRET_1 = np.maximum(0, -LOW)
REGRET_2 = np.maximum(0, HIGH)
assert np.all(LOW <= HIGH)

INK = "#16283C"
BLUE = "#1E6FAF"
PURPLE = "#6947A3"
TEAL = "#0B7A75"
GOLD = "#B47E16"
MUTED = "#526677"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.3,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

fig, axes = plt.subplots(1, 2, figsize=(6.69, 3.78), sharex=True)
fig.patch.set_facecolor("white")
fig.subplots_adjust(left=.09, right=.97, top=.80, bottom=.22, wspace=.30)

for ax in axes:
    ax.set_xlim(.05, .50)
    ax.axvspan(.05, 1 / 6, color="#E8F4F2", zorder=0)
    ax.axvspan(1 / 6, 1 / 3, color="#F2EDF8", zorder=0)
    ax.axvspan(1 / 3, .50, color="#FCEFF0", zorder=0)
    for x in (1 / 6, 1 / 4, 1 / 3):
        ax.axvline(x, color="#AAB7C5", lw=.75, ls=(0, (3, 3)), zorder=1)
    ax.axvline(.20, color=GOLD, lw=1.25, zorder=2)
    ax.set_xticks([.05, 1 / 6, .25, 1 / 3, .4, .5],
                  ["0.05", "1/6", "1/4", "1/3", "0.40", "0.50"])
    ax.tick_params(axis="both", colors=INK, labelsize=7.6, length=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#9AAABC")
    ax.grid(axis="y", color="#E4EAF0", lw=.6, zorder=0)
    ax.set_xlabel("Evaluation threshold  t", color=INK, labelpad=7)

ax = axes[0]
ax.axhline(0, color=INK, lw=.8, zorder=2)
ax.fill_between(T, LOW, HIGH, color="#B9D3E9", alpha=.58, zorder=2)
ax.plot(T, HIGH, color=BLUE, lw=1.8, label="Upper bound", zorder=3)
ax.plot(T, LOW, color=PURPLE, lw=1.8, label="Lower bound", zorder=3)
ax.scatter([.20, .20], [float((1 / 3 - .2) / (3 * .8)),
                         float((1 / 6 - .2) / (3 * .8))],
           color=GOLD, edgecolor="white", lw=.7, s=26, zorder=5)
ax.set_ylim(-.24, .11)
ax.set_ylabel("Net-benefit difference", color=INK)
ax.set_title("A   Compatible fixed-rule contrasts", loc="left", pad=11,
             color=INK, fontsize=9.0, fontweight="bold")
ax.legend(frameon=False, fontsize=7.4, loc="upper right")

ax = axes[1]
ax.plot(T, REGRET_1, color=BLUE, lw=1.8, label="Rule 1", zorder=3)
ax.plot(T, REGRET_2, color=PURPLE, lw=1.8, label="Rule 2", zorder=3)
ax.scatter([.25], [1 / 27], color=TEAL, edgecolor="white", lw=.7, s=32, zorder=5)
ax.annotate("Equal regret\nat t = 1/4", xy=(.25, 1 / 27),
            xytext=(.31, .077), fontsize=7.0, color=INK,
            arrowprops={"arrowstyle": "-", "color": MUTED, "lw": .8})
ax.set_ylim(-.005, .23)
ax.set_ylabel("Worst-case regret", color=INK)
ax.set_title("B   Criterion-dependent choice", loc="left", pad=11,
             color=INK, fontsize=9.0, fontweight="bold")
ax.legend(frameon=False, fontsize=7.4, loc="upper left")

fig.text(.09, .94, "Evidence and preference answer different questions",
         fontsize=10.0, fontweight="bold", color=INK)
fig.text(.09, .89, "The same observationally compatible worlds; two rules held fixed as utility threshold varies",
         fontsize=7.7, color=MUTED)
fig.text(.09, .095,
         "Background: uniform rule 1 · sign ambiguous · uniform rule 2.  Gold: original t = 1/5.",
         fontsize=7.0, color=MUTED)
fig.text(.09, .055,
         "Exact post hoc toy calculation; not a clinical threshold, causal effect, or sampling confidence interval.",
         fontsize=7.0, color=MUTED)

for ext in ("pdf", "svg", "png"):
    fig.savefig(OUT / f"figureS2_threshold_criterion.{ext}",
                dpi=600 if ext == "png" else None)
plt.close(fig)
print("figureS2 saved", OUT)

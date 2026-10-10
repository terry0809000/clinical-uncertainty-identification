"""Figure 3 from verified raw/summary simulation outputs.

Shows Monte Carlo uncertainty in certification and interval width. Direct target
auditing and one-group projection use different non-target evidence; caption
must state that comparison is illustrative, not total-cost matched superiority.
"""

from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1] / "simulation"
OUT = Path(__file__).resolve().parent
main = json.loads((ROOT / "summary_main_10000.json").read_text(encoding="utf-8"))
supp = json.loads((ROOT / "supplement_v2_summary_main_10000.json").read_text(encoding="utf-8"))
assert main["n_repetitions_per_cell"] == supp["reps_per_cell"] == 10000

worlds = (("A_positive", "World A   true ΔNB = +0.0556"),
          ("B_negative", "World B   true ΔNB = −0.0139"))
proj = {(c["world"], c["n_audit"]): c for c in main["cells"] if c["world"] in dict(worlds)}
direct = {(c["world"], c["n_audit"]): c for c in supp["cells"]
          if c["method"] == "direct_two_group_audit" and c["world"] in dict(worlds)}

widths = defaultdict(list)
for file, label, width_name in (("raw_replicates_main_10000.csv", "projection", "interval_width"),
                                ("supplement_v2_raw_main_10000.csv", "direct_two_group_audit", "width")):
    with (ROOT / file).open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["world"] not in dict(worlds):
                continue
            if label == "direct_two_group_audit" and row["method"] != label:
                continue
            value = float(row[width_name])
            if np.isfinite(value):
                widths[(label, row["world"], int(row["n_audit"]))].append(value)


def exact_prop_interval(k: int, n: int, alpha=0.05):
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


COL = {"navy": "#13253D", "projection": "#176BB0", "direct": "#B83B4D",
       "muted": "#5C6B7B", "line": "#DEE6ED"}
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7.5, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
fig, axes = plt.subplots(2, 2, figsize=(6.69, 4.45), sharex="col")
fig.patch.set_facecolor("white")
fig.subplots_adjust(left=0.105, right=0.965, bottom=0.19, top=0.84,
                    wspace=0.26, hspace=0.36)
pos = {0: 0, 50: 1, 200: 2, 800: 3}
handles = []

for col, (world, title) in enumerate(worlds):
    axc, axw = axes[0, col], axes[1, col]
    axc.set_title(("A   " if col == 0 else "B   ") + title, loc="left",
                  fontsize=8.5, fontweight="bold", color=COL["navy"], pad=7)
    axw.set_title(("C   World A · interval width" if col == 0 else
                   "D   World B · interval width"), loc="left",
                  fontsize=8.2, fontweight="bold", color=COL["navy"], pad=7)
    for method, label, cells, color, marker, style in (
        ("projection", "One-group projection", proj, COL["projection"], "o", "-"),
        ("direct_two_group_audit", "Direct two-group audit", direct, COL["direct"], "s", "--"),
    ):
        ns = [n for n in (0, 50, 200, 800) if (world, n) in cells]
        xs = [pos[n] for n in ns]
        cs = [cells[(world, n)] for n in ns]
        vals = [c["metrics"]["correct_sign"]["proportion"] for c in cs]
        ci = [exact_prop_interval(c["metrics"]["correct_sign"]["count"], 10000) for c in cs]
        lo = [100 * (v - a) for v, (a, b) in zip(vals, ci)]
        hi = [100 * (b - v) for v, (a, b) in zip(vals, ci)]
        err = axc.errorbar(xs, np.array(vals) * 100, yerr=np.array([lo, hi]),
                           color=color, marker=marker, linestyle=style,
                           linewidth=1.6, markersize=4.8, capsize=2.3,
                           markeredgecolor="white", markeredgewidth=0.6,
                           label=label, zorder=3)
        if col == 0:
            handles.append(err)
        avgs = [float(np.mean(widths[(method, world, n)])) for n in ns]
        ses = [float(np.std(widths[(method, world, n)], ddof=1) /
                     np.sqrt(len(widths[(method, world, n)]))) for n in ns]
        axw.errorbar(xs, avgs, yerr=np.array(ses) * 1.96,
                     color=color, marker=marker, linestyle=style,
                     linewidth=1.6, markersize=4.8, capsize=2.3,
                     markeredgecolor="white", markeredgewidth=0.6, zorder=3)
    for ax in (axc, axw):
        ax.set_xlim(-0.2, 3.2)
        ax.set_xticks([0, 1, 2, 3], ["0", "50", "200", "800"])
        ax.grid(axis="y", color=COL["line"], lw=0.55)
        ax.set_axisbelow(True)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(COL["line"])
        ax.tick_params(axis="both", labelcolor=COL["muted"], length=0)
    axc.set_ylim(0, 35 if col == 0 else 1.2)
    axw.set_ylim(0, 0.60)
    axw.set_xlabel("Target reference audits, total n", color=COL["navy"])

axes[0, 0].set_ylabel("Correct sign certified (%)", color=COL["navy"])
axes[1, 0].set_ylabel("Mean envelope width", color=COL["navy"])
# Article-like figure title is kept in the manuscript legend, not the artwork.
fig.legend(handles, ["One-group projection", "Direct two-group audit"],
           loc="upper left", bbox_to_anchor=(0.10, 0.965), ncol=2,
           frameon=False, fontsize=7.4, handlelength=2.8)
fig.text(0.055, 0.055,
         "10,000 repetitions per cell. Bars: exact 95% binomial Monte Carlo intervals; width bars: ±1.96 MC SE.\n"
         "Upper-panel vertical scales differ. Projection also uses source paired audits and target proxy cohorts; total costs differ.",
         color=COL["muted"], ha="left", va="bottom", fontsize=7.0)

for ext in ("pdf", "svg", "png"):
    fig.savefig(OUT / f"figure3_simulation.{ext}", dpi=600 if ext == "png" else None)
plt.close(fig)
print("figure3 saved", OUT)

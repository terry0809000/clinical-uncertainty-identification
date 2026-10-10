"""Supplementary exact-population sensitivity surface for clinical uncertainty.

Cells come from the locked V5 rational analytic supplement, not patient data.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / "analytic" / "sensitivity_surface_v5_exact.json"
data = json.loads(SOURCE.read_text(encoding="utf-8"))
assert data["cells"] == 56
dp = [x["decimal"] for x in data["dp_grid"]]
dc = [x["decimal"] for x in data["dc_grid"]]
index = {(r["dp"]["decimal"], r["dc"]["decimal"]): r for r in data["rows"]}
width = np.array([[index[(p, c)]["width"]["decimal"] for c in dc] for p in dp])
status = [[index[(p, c)]["sign_status"] for c in dc] for p in dp]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
fig, ax = plt.subplots(figsize=(6.69, 5.05))
fig.patch.set_facecolor("white")
fig.subplots_adjust(left=0.155, right=0.81, bottom=0.27, top=0.85)
im = ax.imshow(width, cmap="Blues", vmin=0, vmax=0.07,
               origin="lower", interpolation="nearest", aspect="auto")
ax.set_xticks(range(len(dc)), [f"{v:g}" for v in dc])
ax.set_yticks(range(len(dp)), [f"{v:g}" for v in dp])
ax.tick_params(length=0, pad=6, labelsize=8.1)
ax.set_xlabel(r"Channel drift cap  ·  $|Se_1-Se_0|, |Sp_1-Sp_0|$", fontsize=8.8, labelpad=9)
ax.set_ylabel(r"Target risk drift cap  ·  $|p_{1g}-p_{0g}|$", fontsize=8.8, labelpad=10)

for row in range(len(dp)):
    for col in range(len(dc)):
        value = width[row, col]
        mark = "+" if status[row][col] == "strict_positive" else "0"
        ink = "white" if value >= 0.038 else "#16283C"
        ax.text(col, row, mark, ha="center", va="center", color=ink,
                fontsize=11.3, fontweight="bold")
        ax.add_patch(patches.Rectangle((col - .5, row - .5), 1, 1,
                                       facecolor="none", edgecolor="#FFFFFF",
                                       linewidth=.8))
        if mark == "0":
            ax.add_patch(patches.Circle((col, row), radius=.29, fill=False,
                                        edgecolor="white" if value >= .038 else "#A23B47",
                                        linewidth=1.0))

# Outline the original Box 1 assumption pair (dp=.10, dc=.05).
box_row, box_col = dp.index(.10), dc.index(.05)
ax.add_patch(patches.Rectangle((box_col - .47, box_row - .47), .94, .94,
                               fill=False, edgecolor="#D09C2D", linewidth=2.4))
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xlim(-.5, len(dc) - .5)
ax.set_ylim(-.5, len(dp) - .5)

cax = fig.add_axes([.845, .34, .025, .45])
cb = fig.colorbar(im, cax=cax, ticks=[0, .02, .04, .06])
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=7.5, length=0, pad=5)
fig.text(.84, .82, "Set width\n(Delta NB units)", ha="left", va="bottom",
         fontsize=8.0, color="#16283C")

fig.text(.155, .942,
         "Same observed target proxy rates; two sensitivity assumptions vary",
         fontsize=9.0, fontweight="bold", color="#16283C", ha="left")
fig.text(.155, .905,
         "Exact population bounds across 56 locked assumption pairs",
         fontsize=7.8, color="#596B7B", ha="left")
fig.text(.155, .12,
         "+  strictly positive rule contrast     0  zero compatible (outlined)     gold  original Box 1",
         fontsize=7.5, color="#33465A", ha="left")
fig.text(.155, .07,
         "Widths are conditional on stipulated drift caps; these are not confidence intervals or source shares.",
         fontsize=7.3, color="#596B7B", ha="left")

for ext in ("pdf", "svg", "png"):
    fig.savefig(OUT / f"figureS1_sensitivity_surface.{ext}",
                dpi=600 if ext == "png" else None)
plt.close(fig)
print("figureS1 saved", OUT)

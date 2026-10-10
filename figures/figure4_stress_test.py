"""Version-3 synthetic stress test: assumption failures and interaction.

Four heatmap panels show raw-denominator operating characteristics. Only the
upper-left design cell within each panel retains the stipulated audit model.
This is an illustrative secondary analysis, not a clinical frequency estimate.
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import Normalize
import numpy as np

parser = argparse.ArgumentParser(description="Render the exact synthetic stress-test grid as vector heatmap cells.")
parser.add_argument("--data", type=Path,
                    default=Path(__file__).resolve().parents[1] / "legacy" / "simulation" /
                    "supplement_v4_summary_main_10000.json")
parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
args = parser.parse_args()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)
data = json.loads(args.data.read_text(encoding="utf-8"))
assert data["repetitions_per_cell"] == 10000
assert len(data["cells"]) == 8

index = {}
for cell in data["cells"]:
    key = (cell["world"], cell["lambda_proxy_positive"], cell["reference_se"])
    assert key not in index
    index[key] = cell

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.6,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "AIM-Figure4-stress-test-v1",
})
fig, axes = plt.subplots(2, 2, figsize=(6.69, 5.15))
fig.patch.set_facecolor("white")
fig.subplots_adjust(left=0.095, right=0.96, bottom=0.17, top=0.84,
                    wspace=0.38, hspace=0.50)

worlds = [
    ("B_negative", "World B  ·  true ΔNB = −0.0139"),
    ("C_interior_tie", "World C  ·  true ΔNB = 0"),
]
metrics = [
    ("true_covered", "True contrast covered (%)", "Blues", "#174D78"),
    ("wrong_sign", "Incorrect nonzero sign (%)", "Reds", "#9B3545"),
]

for row, (metric, row_label, cmap, accent) in enumerate(metrics):
    for col, (world, world_label) in enumerate(worlds):
        ax = axes[row, col]
        values = np.array([
            [index[(world, lam, se)]["metrics"][metric]["proportion"] * 100
             for se in (1.0, 0.8)]
            for lam in (1.0, 4.0)
        ])
        # Same 0–100% color mapping and imshow tile extent; PDF/SVG vector paths.
        palette = plt.get_cmap(cmap)
        normalize = Normalize(vmin=0, vmax=100)
        for i in range(2):
            for j in range(2):
                ax.add_patch(patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor=palette(normalize(values[i, j])),
                    edgecolor="none", linewidth=0, antialiased=False))
        ax.set_xticks([0, 1], ["Exact", "Imperfect"])
        ax.set_yticks([0, 1], ["Random", "Proxy-selected"])
        ax.tick_params(axis="x", labelsize=7.1, length=0, pad=5)
        ax.tick_params(axis="y", labelsize=7.1, length=0, pad=5)
        if row == 0:
            ax.set_title(("A   " if col == 0 else "B   ") + world_label + "\nTrue contrast covered",
                         loc="left", color="#13253D", fontweight="bold",
                         fontsize=7.9, pad=9)
        else:
            ax.set_title(("C   " if col == 0 else "D   ") + world_label + "\nIncorrect nonzero sign",
                         loc="left", color=accent, fontweight="bold",
                         fontsize=7.9, pad=9)
        for i in range(2):
            for j in range(2):
                value = values[i, j]
                ink = "white" if value >= 55 else "#13253D"
                ax.text(j, i, f"{value:.2f}%", ha="center", va="center",
                        fontsize=10.6, fontweight="bold", color=ink)
        # Distinguish the sole correctly specified control by an explicit border.
        ax.add_patch(patches.Rectangle((-0.48, -0.48), 0.96, 0.96,
                                       fill=False, ec="#182536", lw=1.5))
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(1.5, -0.5)
        ax.set_aspect("equal")

fig.text(0.095, 0.95,
         "Fixed audit n = 800  ·  10,000 repetitions per cell  ·  dark outline: correctly specified audit",
         fontsize=7.5, color="#5C6B7B", va="center")
fig.text(0.095, 0.07,
         "Imperfect reference: sensitivity 0.80, specificity 0.95. Proxy-selected audit: positive proxy weight 4×.\n"
         "The analyst ignores both violations; blue shows coverage, red shows wrong sign; both scales are 0–100%.",
         fontsize=7.0, color="#5C6B7B", va="bottom")

for extension in ("pdf", "svg", "png"):
    metadata = ({"CreationDate": datetime(2026, 9, 25, tzinfo=timezone.utc)}
                if extension == "pdf" else {"Date": "2026-09-25"}
                if extension == "svg" else None)
    fig.savefig(OUT / f"figure4_stress_test.{extension}",
                dpi=600 if extension == "png" else None, metadata=metadata)
plt.close(fig)
print("figure4 saved", OUT)

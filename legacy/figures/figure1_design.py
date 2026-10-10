"""Figure 1: verifiable, editable vector schematic of the proposed evidence design.

All shapes and labels derive from the project's pre-specified analytic/data design.
The interval bars are explicitly schematic; no patient results are plotted.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = Path(__file__).resolve().parent
COL = {
    "navy": "#13253D", "blue": "#1676B8", "cyan": "#00A6C8",
    "teal": "#08756F", "violet": "#6947A8", "amber": "#A67212",
    "coral": "#B83B4D", "ink": "#24364B", "muted": "#5C6B7B",
    "paper": "#F7FAFD", "line": "#D8E2EB", "white": "#FFFFFF",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": COL["paper"],
})

fig, ax = plt.subplots(figsize=(6.69, 5.42))  # 170 mm journal full width
fig.patch.set_facecolor(COL["paper"])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_axis_off()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)


def box(x, y, w, h, face=COL["white"], edge=COL["line"], radius=0.025, lw=0.8):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0.004,rounding_size={radius}",
                       linewidth=lw, edgecolor=edge, facecolor=face,
                       transform=ax.transAxes)
    ax.add_patch(p)
    return p


def txt(x, y, s, size=7.5, color="ink", weight="normal", ha="left", va="center", **kw):
    ax.text(x, y, s, fontsize=size, color=COL.get(color, color),
            fontweight=weight, ha=ha, va=va, transform=ax.transAxes, **kw)


def arrow(x1, y1, x2, y2, color="muted", lw=1.2, style="-|>", **kw):
    p = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                        mutation_scale=8, linewidth=lw,
                        color=COL.get(color, color), transform=ax.transAxes, **kw)
    ax.add_patch(p)


# Dark header is a four-stage process graphic; the figure title is in the manuscript legend.
box(0.018, 0.842, 0.964, 0.139, face=COL["navy"], edge=COL["navy"], radius=0.026)
stages = [
    (0.047, "TARGET", "Who, when and under which care?", "blue"),
    (0.292, "EVIDENCE", "What is observed or audited?", "cyan"),
    (0.542, "IDENTIFICATION", "What can the evidence imply?", "violet"),
    (0.817, "ACTION", "Which rule is justified?", "teal"),
]
for x, title, sub, key in stages:
    txt(x, 0.927, title, 9.3, "white", "bold")
    txt(x, 0.886, sub, 6.8, "#DCEAF3")
    ax.add_patch(Rectangle((x, 0.86), 0.085, 0.007,
                           facecolor=COL[key], edgecolor="none", transform=ax.transAxes))
for x in (0.255, 0.508, 0.785):
    txt(x, 0.925, "→", 9.5, "#DCEAF3")

# Three main panels.
panels = [(0.023, 0.275, 0.302, 0.54), (0.348, 0.275, 0.302, 0.54),
          (0.674, 0.275, 0.302, 0.54)]
for x, y, w, h in panels:
    box(x, y, w, h)

txt(0.046, 0.784, "01  DEFINE TARGET", 8.1, "blue", "bold")
txt(0.371, 0.784, "02  NAME UNCERTAINTY", 8.1, "violet", "bold")
txt(0.697, 0.784, "03  EVALUATE ACTION", 8.1, "teal", "bold")

# Panel 1: sites, patient timeline, reference and proxy outcome.
box(0.05, 0.675, 0.105, 0.055, face="#E6F3FA", edge="#B6D9E8", radius=0.011)
box(0.195, 0.675, 0.105, 0.055, face="#E7F6F3", edge="#B4DED8", radius=0.011)
txt(0.102, 0.702, "SOURCE", 7.2, "blue", "bold", "center")
txt(0.247, 0.702, "TARGET", 7.2, "teal", "bold", "center")
arrow(0.158, 0.703, 0.19, 0.703, "muted", 1.0)
txt(0.174, 0.748, "shift", 7.0, "muted", ha="center")
txt(0.102, 0.649, "paired audit", 7.0, "blue", ha="center")
txt(0.247, 0.649, "proxy cohort", 7.0, "teal", ha="center")

ax.add_patch(Circle((0.082, 0.566), 0.014, facecolor=COL["blue"], transform=ax.transAxes))
box(0.067, 0.514, 0.03, 0.037, face=COL["blue"], edge=COL["blue"], radius=0.011)
txt(0.108, 0.566, "eligible patient", 7.1)
ax.plot([0.07, 0.282], [0.481, 0.481], color=COL["line"], lw=1.6,
        transform=ax.transAxes, clip_on=False)
for x, k in ((0.092, "blue"), (0.279, "coral")):
    ax.add_patch(Circle((x, 0.481), 0.008, facecolor=COL[k],
                        edgecolor=COL["white"], lw=0.5, transform=ax.transAxes))
txt(0.071, 0.451, "index + score", 7.0, "blue")
txt(0.299, 0.451, "endpoint", 7.0, "coral", ha="right")

box(0.048, 0.354, 0.118, 0.055, face="#FFF4DB", edge="#EBD5A3", radius=0.012)
box(0.183, 0.354, 0.119, 0.055, face="#FDECEE", edge="#EEC5CD", radius=0.012)
txt(0.107, 0.381, "routine proxy", 7.0, "#8C6115", "bold", "center")
txt(0.243, 0.381, "audit reference", 7.0, "coral", "bold", "center")
txt(0.049, 0.321, "Reference needs an audit.", 7.0, "muted")
txt(0.049, 0.302, "Care policy fixes the target.", 7.0, "muted")

# Panel 2: four scientifically distinct questions, each with redundant symbol encoding.
layers = [
    ("blue", "o", "Future outcome", "Variation after risk is fixed"),
    ("amber", "D", "Data process", "Measurement and labels"),
    ("violet", "s", "Model estimate", "Fitting and specification"),
    ("teal", "^", "Population / policy", "Setting and care pathway"),
]
for i, (key, marker, label, sub) in enumerate(layers):
    yy = 0.68 - i * 0.104
    box(0.371, yy - 0.041, 0.254, 0.087, face="#FBFCFE", edge=COL["line"], radius=0.011)
    ax.plot([0.395], [yy + 0.012], marker=marker, markersize=6.2,
            markerfacecolor=COL[key], markeredgecolor=COL[key],
            linestyle="none", transform=ax.transAxes)
    txt(0.417, yy + 0.019, label, 7.3, "ink", "bold")
    txt(0.417, yy - 0.018, sub, 7.0, "muted")
txt(0.371, 0.294, "Mechanisms may overlap.", 7.0, "muted")

# Panel 3: locked rules, disagreement, feasible interval, and audit.
txt(0.697, 0.749, "Which rule is preferable at t?", 7.0, "muted")
box(0.698, 0.676, 0.115, 0.054, face="#E6F3FA", edge="#B6D9E8", radius=0.011)
box(0.839, 0.676, 0.11, 0.054, face="#EEEAF8", edge="#CDC2E9", radius=0.011)
txt(0.755, 0.703, "Rule 1", 7.0, "blue", "bold", "center")
txt(0.894, 0.703, "Rule 2", 7.0, "violet", "bold", "center")
arrow(0.777, 0.669, 0.813, 0.617, "blue", 1.0)
arrow(0.88, 0.669, 0.842, 0.617, "violet", 1.0)
box(0.724, 0.56, 0.218, 0.064, face="#F4F2FA", edge="#D9D2EB", radius=0.012)
txt(0.833, 0.595, "Rule disagreement", 7.5, "ink", "bold", "center")
txt(0.833, 0.571, "ranking depends on reference", 7.0, "muted", ha="center")

txt(0.698, 0.521, "Net-benefit difference", 7.2, "ink", "bold")
ax.plot([0.705, 0.949], [0.457, 0.457], color=COL["line"], lw=1.1, transform=ax.transAxes)
ax.plot([0.826, 0.826], [0.439, 0.475], color=COL["ink"], lw=0.8, transform=ax.transAxes)
txt(0.826, 0.424, "0", 7.0, "muted", ha="center")
ax.plot([0.735, 0.916], [0.457, 0.457], color="#C7DBEC", lw=9.2,
        solid_capstyle="round", transform=ax.transAxes)
ax.plot([0.751, 0.893], [0.457, 0.457], color=COL["coral"], lw=5.2,
        solid_capstyle="round", transform=ax.transAxes)
txt(0.828, 0.486, "identified set", 7.0, "coral", ha="center")
txt(0.828, 0.399, "pale = sampling envelope", 7.0, "muted", ha="center")
box(0.703, 0.31, 0.246, 0.062, face="#E7F6F3", edge="#B4DED8", radius=0.012)
txt(0.826, 0.343, "Sampled G+ audit", 7.2, "teal", "bold", "center")
txt(0.826, 0.320, "may narrow or flag model conflict", 7.0, "muted", ha="center")

# Footer: a concise experimental programme, clearly labelled as design.
box(0.023, 0.041, 0.953, 0.205, face="#ECF2F7", edge="#DAE5EE", radius=0.02)
txt(0.045, 0.22, "EVIDENCE PROGRAM  ·  review / analytic / simulation", 8.2, "navy", "bold")
footer = [
    (0.047, "01", "Method evidence map", "Target, assumptions, coverage", "blue"),
    (0.355, "02", "Paired worlds", "Same data; reversed rank", "violet"),
    (0.666, "03", "Locked simulation", "Coverage and failure modes", "teal"),
]
for x, num, title, sub, key in footer:
    box(x, 0.073, 0.286, 0.11, face=COL["white"], edge=COL["line"], radius=0.014)
    ax.add_patch(Circle((x + 0.029, 0.128), 0.019, facecolor=COL[key], transform=ax.transAxes))
    txt(x + 0.029, 0.129, num, 7.0, "white", "bold", "center")
    txt(x + 0.057, 0.143, title, 7.25, "ink", "bold")
    txt(x + 0.057, 0.107, sub, 7.0, "muted")
txt(0.975, 0.018, "Conceptual design · interval schematic", 7.0, "muted", ha="right")

for ext in ("pdf", "svg", "png"):
    fig.savefig(OUT / f"figure1_design.{ext}", dpi=600 if ext == "png" else None)
plt.close(fig)
print("exported", *(str(OUT / f"figure1_design.{e}") for e in ("pdf", "svg", "png")))

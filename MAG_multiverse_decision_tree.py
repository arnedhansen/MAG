"""
AOC Multiverse — Traditional vs Multiverse Analysis Comparison
Single presentation figure showing one arbitrary analysis path (traditional)
alongside all specifications (multiverse) using trunk-and-reconverge layout.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_GH = ("/Users/Arne/Documents/GitHub/MAG/multiverse/decision_tree/"
          "MAG_multiverse_decision_tree_traditional_vs_multiverse.png")
OUT_SC = ("/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/"
          "figures/multiverse/decision_tree/"
          "MAG_multiverse_decision_tree_traditional_vs_multiverse.png")

# ── Palette (minimal, near-monochrome) ───────────────────────────────────
CHARCOAL  = "#374151"
SLATE     = "#64748B"
LIGHT_BG  = "#F1F5F9"
BORDER    = "#94A3B8"
GUIDE     = "#E2E8F0"

# ── 7 decisions (display order) ──────────────────────────────────────────
decisions = [
    {"label": "Latency\nWindow",               "trad": "1000–2000 ms",   "multi": ["0–1000 ms", "1000–2000 ms", "0–2000 ms"]},
    {"label": "Electrodes",          "trad": "Occipital",   "multi": ["Posterior", "Occipital"]},
    {"label": "Spectral\nParameterization",    "trad": "No SpecParam", "multi": ["SpecParam", "No SpecParam"]},
    {"label": "EEG\nBaseline",       "trad": "Raw",         "multi": ["Raw", "%"]},
    {"label": "Alpha\nBand",         "trad": "IAF",         "multi": ["Canonical", "IAF"]},
    {"label": "Gaze\nMeasure",       "trad": "Gaze Deviation",   "multi": ["Gaze Deviation", "SPL", "Velocity", "BCEA"]},
    {"label": "Gaze\nBaseline",      "trad": "Raw",         "multi": ["Raw", "%"]},
]
N = len(decisions)

# ── Figure & axes ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(26, 19), facecolor="white")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# ── Spatial zones ────────────────────────────────────────────────────────
LABEL_CX    = 0.055
TRAD_CX     = 0.195
DIVIDER_X   = 0.30
MULTI_L     = 0.37
MULTI_R     = 0.97
MULTI_TRUNK = (MULTI_L + MULTI_R) / 2

Y_ROOT   = 0.925
Y_TOP    = 0.845
Y_BOT    = 0.155
Y_LEAF   = 0.055
NODE_H   = 0.048
GAP      = 0.012
BOX_PAD  = 0.004

y_levels = np.linspace(Y_TOP, Y_BOT, N)


# ── Drawing helpers ──────────────────────────────────────────────────────
def rbox(cx, cy, w, h, fc, text, fs=10, tc="white", bold=False,
         ec=None, lw=2):
    ax.add_patch(mpatches.FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=mpatches.BoxStyle.Round(pad=BOX_PAD),
        facecolor=fc, edgecolor=ec or fc, linewidth=lw,
        zorder=3, transform=ax.transAxes, clip_on=False))
    if text:
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
                color=tc, fontweight="bold" if bold else "normal",
                fontfamily="sans-serif", transform=ax.transAxes,
                zorder=4, linespacing=1.15)


def scurve(x1, y1, x2, y2, col=SLATE, lw=1.6, alpha=0.35):
    t = np.linspace(0, 1, 40)
    s = 3 * t**2 - 2 * t**3
    ax.plot(x1 + (x2 - x1) * s, y1 + (y2 - y1) * t,
            color=col, lw=lw, alpha=alpha, solid_capstyle="round",
            transform=ax.transAxes, zorder=1)


def varrow(x, y1, y2, col=SLATE, lw=1.6):
    ax.annotate(
        "",
        xy=(x, y2),
        xytext=(x, y1),
        xycoords=ax.transAxes,
        textcoords=ax.transAxes,
        arrowprops=dict(
            arrowstyle="-|>",
            color=col,
            lw=lw,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=18,
        ),
        zorder=5,
    )


# ── Titles ───────────────────────────────────────────────────────────────
ax.text(TRAD_CX, 0.98, "Traditional Analysis",
        ha="center", va="center", fontsize=30, fontweight="bold",
        color=CHARCOAL, fontfamily="sans-serif")
ax.text(MULTI_TRUNK, 0.98, "Multiverse Analysis",
        ha="center", va="center", fontsize=30, fontweight="bold",
        color=CHARCOAL, fontfamily="sans-serif")

# ── Divider ──────────────────────────────────────────────────────────────
ax.plot([DIVIDER_X, DIVIDER_X], [0.015, 0.965],
        color=GUIDE, lw=2, transform=ax.transAxes, zorder=0)

# ── Start nodes (wide rounded rectangle, not a circle) ──────────────────
start_w, start_h = 0.07, 0.030
rbox(TRAD_CX, Y_ROOT, start_w, start_h, CHARCOAL, "Start",
     fs=12, tc="white", bold=True)
rbox(MULTI_TRUNK, Y_ROOT, start_w, start_h, CHARCOAL, "Start",
     fs=12, tc="white", bold=True)

# ── Multiverse trunk (background dashed line) ───────────────────────────
ax.plot([MULTI_TRUNK, MULTI_TRUNK],
        [Y_ROOT - start_h / 2 - 0.005, y_levels[-1] - NODE_H / 2 - GAP - 0.005],
        color=GUIDE, lw=2.5, ls="--", alpha=0.55,
        transform=ax.transAxes, zorder=0)

# ── Draw each decision level ────────────────────────────────────────────
running = 1

for i, dec in enumerate(decisions):
    y = y_levels[i]
    opts = dec["multi"]
    n_opts = len(opts)
    running *= n_opts

    # ── Decision label (far left) ──
    ax.text(LABEL_CX, y, dec["label"], ha="center", va="center",
            fontsize=18, color=CHARCOAL, fontweight="bold",
            fontfamily="sans-serif", transform=ax.transAxes,
            linespacing=1.1)

    # Dotted guideline across
    ax.plot([LABEL_CX + 0.06, TRAD_CX - 0.05], [y, y],
            color=GUIDE, lw=1, ls=":", transform=ax.transAxes, zorder=0)

    # ── TRADITIONAL: single box ──
    prev_bot = ((Y_ROOT - start_h / 2 - BOX_PAD)
                if i == 0 else (y_levels[i - 1] - NODE_H / 2 - BOX_PAD))
    rbox(TRAD_CX, y, 0.085, NODE_H, LIGHT_BG, dec["trad"],
         fs=13.5, tc=CHARCOAL, ec=BORDER, lw=1.8)
    varrow(TRAD_CX, prev_bot, y + NODE_H / 2 + BOX_PAD)

    # ── MULTIVERSE: trunk → fan out → boxes → fan in ──
    span = MULTI_R - MULTI_L
    if n_opts == 2:
        half = span * 0.25
        opt_xs = [MULTI_TRUNK - half, MULTI_TRUNK + half]
    elif n_opts <= 4:
        m = span * 0.06
        opt_xs = list(np.linspace(MULTI_L + m, MULTI_R - m, n_opts))
    else:
        m = span * 0.03
        opt_xs = list(np.linspace(MULTI_L + m, MULTI_R - m, n_opts))

    bw = 0.10 if n_opts <= 2 else (0.095 if n_opts <= 4 else 0.08)

    split_y = y + NODE_H / 2 + GAP + 0.004
    conv_y  = y - NODE_H / 2 - GAP - 0.004

    # Trunk segment from previous converge to this split
    if i == 0:
        ax.plot([MULTI_TRUNK, MULTI_TRUNK],
                [Y_ROOT - start_h / 2, split_y],
                color=SLATE, lw=1.5, alpha=0.35, transform=ax.transAxes, zorder=1)
    else:
        prev_conv = y_levels[i - 1] - NODE_H / 2 - GAP - 0.004
        ax.plot([MULTI_TRUNK, MULTI_TRUNK], [prev_conv, split_y],
                color=SLATE, lw=1.5, alpha=0.35, transform=ax.transAxes, zorder=1)

    # Split dot
    ax.plot(MULTI_TRUNK, split_y, "o", color=CHARCOAL, markersize=5,
            zorder=5, transform=ax.transAxes)

    for ox, label in zip(opt_xs, opts):
        scurve(MULTI_TRUNK, split_y, ox, y + NODE_H / 2)
        scurve(ox, y - NODE_H / 2, MULTI_TRUNK, conv_y)
        rbox(ox, y, bw, NODE_H, LIGHT_BG, label,
             fs=13 if n_opts <= 4 else 12, tc=CHARCOAL, ec=BORDER, lw=1.8)

    # Converge dot
    ax.plot(MULTI_TRUNK, conv_y, "o", color=CHARCOAL, markersize=5,
            zorder=5, transform=ax.transAxes)

    # Universe count (right margin, outside plot area)
    ax.text(1.02, y + 0.007, f"×{n_opts}", ha="left", va="center",
            fontsize=15, color=CHARCOAL, fontweight="bold",
            fontfamily="sans-serif", transform=ax.transAxes)
    ax.text(1.02, y - 0.017, f"= {running:,}", ha="left", va="center",
            fontsize=14, color=SLATE,
            fontfamily="sans-serif", transform=ax.transAxes)

# ── Leaf level ───────────────────────────────────────────────────────────
# Label
ax.text(LABEL_CX, Y_LEAF, "Reported\nAnalysis", ha="center", va="center",
        fontsize=15, color=CHARCOAL, fontweight="bold",
        fontfamily="sans-serif", transform=ax.transAxes, linespacing=1.1)
ax.plot([LABEL_CX + 0.06, TRAD_CX - 0.05], [Y_LEAF, Y_LEAF],
        color=GUIDE, lw=1, ls=":", transform=ax.transAxes, zorder=0)

# Traditional: arrow → one tiny box
rbox(TRAD_CX, Y_LEAF, 0.022, 0.016, LIGHT_BG, "", ec=BORDER, lw=1.2)
varrow(TRAD_CX, y_levels[-1] - NODE_H / 2 - BOX_PAD, Y_LEAF + 0.008 + BOX_PAD)
ax.text(TRAD_CX, Y_LEAF - 0.025, "1 analysis", ha="center", va="center",
        fontsize=14, color=SLATE, fontweight="bold",
        fontfamily="sans-serif", transform=ax.transAxes)

# Multiverse: fan out to many tiny boxes
last_conv = y_levels[-1] - NODE_H / 2 - GAP - 0.004
n_tiny = 24
tiny_xs = np.linspace(MULTI_L + 0.01, MULTI_R - 0.01, n_tiny)
for tx in tiny_xs:
    scurve(MULTI_TRUNK, last_conv, tx, Y_LEAF + 0.010,
           col=SLATE, lw=0.7, alpha=0.18)
    rbox(tx, Y_LEAF, 0.016, 0.013, LIGHT_BG, "", ec=BORDER, lw=0.7)

ax.text(MULTI_TRUNK, Y_LEAF - 0.025, f"{running:,} analyses per task",
        ha="center", va="center", fontsize=15, color=CHARCOAL,
        fontweight="bold", fontfamily="sans-serif", transform=ax.transAxes)

# ── Save ─────────────────────────────────────────────────────────────────
import os
for path in (OUT_GH, OUT_SC):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {path}")
plt.close(fig)

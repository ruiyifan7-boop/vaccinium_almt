#!/usr/bin/env python3
"""Figure 1 v2 (taller, no right bracket): tree | name | eco | fam | bars
Eco grouping shown via:
  - colored Eco strip per species (3rd column)
  - pastel background band behind each bar (per ecology)
  - bottom legend listing the 4 ecology categories
"""
import os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

proj = os.path.expanduser("~/projects/vaccinium_almt")
fig_dir = f"{proj}/manuscript/figures"
os.makedirs(fig_dir, exist_ok=True)

tree_def = (
    ("Vvin", ("Ath", "Mdo")),
    (("Csin", "Aer"),
     ("Rsim",
      (("Vmac", "Vmic"),
       ("Vcor", ("Vdar", "Vduc")))))
)

species_info = {
    "Vvin": ("V. vinifera",    12, "Vitaceae",      "Moderate-pH",        ""),
    "Ath":  ("A. thaliana",    16, "Brassicaceae",  "Neutral",            ""),
    "Mdo":  ("M. domestica",   26, "Rosaceae",      "Moderate-pH",        ""),
    "Csin": ("C. sinensis",    17, "Theaceae",      "Acid-tolerant",      ""),
    "Aer":  ("A. eriantha",    23, "Actinidiaceae", "Acid-tolerant",      ""),
    "Rsim": ("R. simsii",      11, "Ericaceae",     "Acid-soil obligate", ""),
    "Vmac": ("V. macrocarpon", 10, "Ericaceae",     "Acid-soil obligate", ""),
    "Vmic": ("V. microcarpum",  9, "Ericaceae",     "Acid-soil obligate", ""),
    "Vcor": ("V. corymbosum",  13, "Ericaceae",     "Acid-soil obligate", "(4x)"),
    "Vdar": ("V. darrowii",     8, "Ericaceae",     "Acid-soil obligate", ""),
    "Vduc": ("V. duclouxii",    9, "Ericaceae",     "Acid-soil obligate", ""),
}

lit_counts = {"Ath": 14, "Csin": 16, "Mdo": 25}

eco_color = {
    "Acid-soil obligate": "#C0392B",
    "Acid-tolerant":      "#27AE60",
    "Moderate-pH":        "#2980B9",
    "Neutral":            "#E67E22",
}
eco_bg = {
    "Acid-soil obligate": "#FADBD8",
    "Acid-tolerant":      "#D5F5E3",
    "Moderate-pH":        "#D6EAF8",
    "Neutral":            "#FAE5D3",
}
family_color = {
    "Ericaceae":     "#922B21",
    "Theaceae":      "#27AE60",
    "Actinidiaceae": "#16A085",
    "Rosaceae":      "#2980B9",
    "Vitaceae":      "#5DADE2",
    "Brassicaceae":  "#E67E22",
}

tip_y = {}
def collect_tips(node):
    if isinstance(node, str):
        tip_y[node] = len(tip_y)
        return tip_y[node]
    return sum(collect_tips(c) for c in node) / len(node)

def tree_depth(node, d=0):
    if isinstance(node, str): return d
    return max(tree_depth(c, d+1) for c in node)

collect_tips(tree_def)
D_MAX = tree_depth(tree_def)

branch_segments = []
def draw_tree(node, depth=0):
    if isinstance(node, str):
        return tip_y[node]
    children_y = [draw_tree(c, depth+1) for c in node]
    y_min, y_max = min(children_y), max(children_y)
    branch_segments.append((depth, depth, y_min, y_max))
    for c, cy in zip(node, children_y):
        end_x = D_MAX if isinstance(c, str) else depth + 1
        branch_segments.append((depth, end_x, cy, cy))
    return sum(children_y) / len(children_y)
draw_tree(tree_def)

species_order = sorted(tip_y.keys(), key=lambda k: tip_y[k])
BACKSLASH = chr(92)

y_tick_labels = []
for abbr in species_order:
    full_name, _, _, _, special = species_info[abbr]
    safe = full_name.replace('. ', '.' + BACKSLASH + ',')
    lbl = '$' + BACKSLASH + 'it{' + safe + '}$'
    if special:
        lbl = lbl + '  ' + special
    y_tick_labels.append(lbl)

# ============================================================
# Layout: 5 panels, NO right bracket panel; taller aspect ratio
# ============================================================
fig = plt.figure(figsize=(10.5, 7.5))
gs = fig.add_gridspec(
    1, 5,
    width_ratios=[1.2, 2.5, 0.3, 0.3, 5.0],
    wspace=0.06,
)
ax_tree = fig.add_subplot(gs[0])
ax_name = fig.add_subplot(gs[1])
ax_eco  = fig.add_subplot(gs[2])
ax_fam  = fig.add_subplot(gs[3])
ax_bar  = fig.add_subplot(gs[4])

N = len(species_order)
YLIM = (-0.7, N - 0.3)

def style_yaxis(ax):
    ax.set_ylim(YLIM)
    ax.invert_yaxis()
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)
    for s in ('top', 'right', 'left', 'bottom'):
        ax.spines[s].set_visible(False)

# ---------- Panel 1: tree ----------
for x0, x1, y0, y1 in branch_segments:
    ax_tree.plot([x0, x1], [y0, y1], color='#222222', linewidth=1.3, zorder=3)
ax_tree.set_xlim(-0.2, D_MAX + 0.1)
style_yaxis(ax_tree)

# ---------- Panel 2: italic species names ----------
ax_name.set_xlim(0, 1)
for i, abbr in enumerate(species_order):
    y = tip_y[abbr]
    ax_name.text(0.02, y, y_tick_labels[i],
                 ha='left', va='center', fontsize=12)
style_yaxis(ax_name)

# ---------- Panel 3: Eco strip ----------
ax_eco.set_xlim(0, 1)
for abbr in species_order:
    y = tip_y[abbr]
    edaphic = species_info[abbr][3]
    ax_eco.add_patch(Rectangle((0.15, y - 0.35), 0.7, 0.7,
                                facecolor=eco_color[edaphic],
                                edgecolor='black', linewidth=0.5))
# Column header above strip
ax_eco.text(0.5, -0.55, "Eco", ha='center', va='center',
            fontsize=10, fontweight='bold', color='#444444')
style_yaxis(ax_eco)

# ---------- Panel 4: Fam strip ----------
ax_fam.set_xlim(0, 1)
for abbr in species_order:
    y = tip_y[abbr]
    family = species_info[abbr][2]
    ax_fam.add_patch(Rectangle((0.15, y - 0.35), 0.7, 0.7,
                                facecolor=family_color[family],
                                edgecolor='black', linewidth=0.5))
ax_fam.text(0.5, -0.55, "Fam", ha='center', va='center',
            fontsize=10, fontweight='bold', color='#444444')
style_yaxis(ax_fam)

# ---------- Panel 5: bars with ecology background bands ----------
ax_bar.set_xlim(0, 31)
ax_bar.set_ylim(YLIM)
ax_bar.invert_yaxis()

for abbr in species_order:
    y = tip_y[abbr]
    edaphic = species_info[abbr][3]
    ax_bar.add_patch(Rectangle((0, y - 0.45), 31, 0.9,
                                facecolor=eco_bg[edaphic], alpha=0.5,
                                edgecolor='none', zorder=0))

for abbr in species_order:
    y = tip_y[abbr]
    count = species_info[abbr][1]
    family = species_info[abbr][2]
    ax_bar.barh(y, count, height=0.6, color=family_color[family],
                edgecolor='black', linewidth=0.6, zorder=3)
    ax_bar.text(count + 0.4, y, str(count), ha='left', va='center',
                fontsize=11, fontweight='bold', zorder=4)

for abbr, lit_n in lit_counts.items():
    if abbr in tip_y:
        y = tip_y[abbr]
        ax_bar.plot([lit_n, lit_n], [y - 0.3, y + 0.3],
                    color='#555555', linestyle='--', linewidth=1.2,
                    alpha=0.85, zorder=2)
        ax_bar.text(lit_n, y - 0.5, "lit:" + str(lit_n), fontsize=8,
                    color='#555555', ha='center', va='top',
                    style='italic', zorder=4)

ax_bar.set_xlabel("Number of ALMT genes per genome", fontsize=12)
ax_bar.set_xticks(np.arange(0, 31, 5))
ax_bar.tick_params(axis='y', left=False, labelleft=False)
ax_bar.grid(True, axis='x', alpha=0.3, linestyle=':', zorder=1)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.spines['left'].set_visible(False)

# ============================================================
# Bottom legend: two rows (Eco + Family)
# ============================================================
eco_handles = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black',
                          linewidth=0.5) for c in eco_color.values()]
fam_handles = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black',
                          linewidth=0.5) for c in family_color.values()]

leg1 = fig.legend(eco_handles, list(eco_color.keys()),
                  loc='lower center', bbox_to_anchor=(0.5, -0.005),
                  ncol=4, fontsize=10, frameon=True, framealpha=0.9,
                  title="Ecology (Eco strip + background band)",
                  title_fontsize=10)
fig.add_artist(leg1)
fig.legend(fam_handles, list(family_color.keys()),
           loc='lower center', bbox_to_anchor=(0.5, -0.075),
           ncol=6, fontsize=10, frameon=True, framealpha=0.9,
           title="Family (Fam strip + bar colour)",
           title_fontsize=10)

fig.suptitle("Genome-wide ALMT family size across 11 angiosperm species",
             fontsize=13.5, fontweight='bold', y=0.97)

plt.tight_layout(rect=[0, 0.1, 1, 0.95])

for fmt in ("svg", "pdf", "png"):
    out = fig_dir + "/Figure1_ALMT_count." + fmt
    plt.savefig(out, format=fmt, dpi=300 if fmt == "pdf" else 200,
                bbox_inches='tight')
    print("  wrote " + out)
plt.close()

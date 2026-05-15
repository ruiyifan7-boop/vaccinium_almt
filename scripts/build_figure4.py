#!/usr/bin/env python3
"""Figure 4: synteny + duplication mode. Based on v6 final, +SVG output."""
import os, re, random
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.path import Path
import matplotlib.patches as mpatches
import numpy as np
from collections import defaultdict

proj = os.path.expanduser("~/projects/vaccinium_almt")
fig_dir = f"{proj}/manuscript/figures"
os.makedirs(fig_dir, exist_ok=True)

almt_ids = {}
for ln in open(f"{proj}/results/p2_almt/combined/ALL_ALMT_verified_ids.txt"):
    m = re.match(r'^([A-Z][a-z]+)_(.+)$', ln.strip())
    if m: almt_ids[m.group(2)] = m.group(1)

def load_gff(p):
    d = {}
    for line in open(p):
        c = line.strip().split("\t")
        if len(c) < 4: continue
        d[c[1]] = (c[0], int(c[2]), int(c[3]))
    return d

vduc_gff = load_gff(f"{proj}/results/p5_evolution/mcscanx/Vduc_single/Vduc.gff")
mdo_gff  = load_gff(f"{proj}/results/p5_evolution/mcscanx/Mdo_single/Mdo.gff")
vdar_gff = load_gff(f"{proj}/results/p5_evolution/mcscanx_vdar_vduc/Vdar.gff")

VDAR_CHRS = [f"Va_{i}" for i in range(151, 163)]
VDUC_CHRS = [f"Vd{i}" for i in range(1, 13)]
mdo_chr_genes = defaultdict(int)
for gid, (chrom, s, e) in mdo_gff.items():
    if chrom.startswith("Md_"): mdo_chr_genes[chrom] += 1
MDO_CHRS = sorted(mdo_chr_genes.keys(), key=lambda c: -mdo_chr_genes[c])[:17]

def chr_lengths(gff, chr_list):
    lens = {}
    for gid, (chrom, s, e) in gff.items():
        if chrom in chr_list: lens[chrom] = max(lens.get(chrom, 0), e)
    return [(c, lens.get(c, 0)) for c in chr_list]

vduc_lens = chr_lengths(vduc_gff, VDUC_CHRS)
mdo_lens  = chr_lengths(mdo_gff, MDO_CHRS)
vdar_lens = chr_lengths(vdar_gff, VDAR_CHRS)

def parse_coll(p, p1, p2):
    out = []; ct = None
    for line in open(p):
        if line.startswith("## Alignment"):
            m = re.search(r'Alignment .*?(\S+)&(\S+)', line)
            if m:
                c1, c2 = m.group(1), m.group(2)
                t1 = p1 if c1.startswith(p1) else (p2 if c1.startswith(p2) else None)
                t2 = p1 if c2.startswith(p1) else (p2 if c2.startswith(p2) else None)
                ct = "cross" if t1 and t2 and t1 != t2 else None
        elif line.strip() and not line.startswith("#") and ct:
            parts = line.split()
            if len(parts) >= 4: out.append((parts[-3], parts[-2]))
    return out

vduc_mdo  = parse_coll(f"{proj}/results/p5_evolution/mcscanx/Vduc_Mdo/VducMdo.collinearity", "Vd", "Md")
vdar_vduc = parse_coll(f"{proj}/results/p5_evolution/mcscanx_vdar_vduc/VdarVduc.collinearity", "Va", "Vd")

TOTAL_W = 100.0
def make_pos(lens_list, total_w, gap_frac=0.012):
    n = len(lens_list); total = sum(L for _, L in lens_list)
    gap = total_w * gap_frac; avail = total_w - gap * (n - 1)
    pos, cum = {}, 0
    for c, L in lens_list:
        w = avail * L / total
        pos[c] = (cum, cum + w); cum += w + gap
    return pos

vdar_pos = make_pos(vdar_lens, TOTAL_W)
vduc_pos = make_pos(vduc_lens, TOTAL_W)
mdo_pos  = make_pos(mdo_lens, TOTAL_W)

Y_MDO, Y_VDUC, Y_VDAR = 6.0, 3.0, 0.0
TH = 0.45

def gene_pos(gff, cp, lens_dict):
    out = {}
    for gid, (chrom, s, e) in gff.items():
        if chrom in cp:
            xs, xe = cp[chrom]
            cl = next((L for c, L in lens_dict if c == chrom), 0)
            if cl == 0: continue
            mid = (s + e) / 2
            out[gid] = xs + (mid / cl) * (xe - xs)
    return out

vdar_gx = gene_pos(vdar_gff, vdar_pos, vdar_lens)
vduc_gx = gene_pos(vduc_gff, vduc_pos, vduc_lens)
mdo_gx  = gene_pos(mdo_gff,  mdo_pos,  mdo_lens)

fig = plt.figure(figsize=(15, 9.5))
gs = fig.add_gridspec(2, 1, height_ratios=[1.6, 1.0], hspace=0.32)
ax_rib = fig.add_subplot(gs[0])
ax_dup = fig.add_subplot(gs[1])
ax_rib.set_xlim(-4, TOTAL_W + 3); ax_rib.set_ylim(-1.5, 7.5)

track_colors = {"Mdo":("#3498DB","#E8F4FD"), "Vduc":("#1F618D","#D6EAF8"), "Vdar":("#C0392B","#FADBD8")}

def draw_track(ax, cp, y, h, ec, fc, label, lo=-3.5):
    for i, (c, (xs, xe)) in enumerate(cp.items(), 1):
        rect = Rectangle((xs, y - h/2), xe - xs, h,
                         facecolor=fc, edgecolor=ec, linewidth=0.8, zorder=5)
        ax.add_patch(rect)
        ax.text((xs+xe)/2, y, str(i), ha='center', va='center',
                fontsize=10, fontweight='bold', color=ec, zorder=6)
    ax.text(lo, y, label, ha='right', va='center',
            fontsize=13, style='italic', fontweight='bold', color=ec)

draw_track(ax_rib, mdo_pos,  Y_MDO,  TH, *track_colors["Mdo"],  'M. domestica\n(17 chr)')
draw_track(ax_rib, vduc_pos, Y_VDUC, TH, *track_colors["Vduc"], 'V. duclouxii\n(12 chr)')
draw_track(ax_rib, vdar_pos, Y_VDAR, TH, *track_colors["Vdar"], 'V. darrowii\n(12 chr)')

def bezier(ax, x1, y1, x2, y2, color, alpha, lw=0.4, zorder=2):
    midy = (y1 + y2) / 2
    verts = [(x1, y1), (x1, midy), (x2, midy), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    patch = mpatches.PathPatch(Path(verts, codes), edgecolor=color, facecolor='none',
                                linewidth=lw, alpha=alpha, zorder=zorder)
    ax.add_patch(patch)

random.seed(42)
for g1, g2 in random.sample(vduc_mdo, min(2500, len(vduc_mdo))):
    if g1 in vduc_gx and g2 in mdo_gx:
        bezier(ax_rib, vduc_gx[g1], Y_VDUC + TH/2, mdo_gx[g2], Y_MDO - TH/2, '#A5B5C0', 0.13, 0.25, 2)
    elif g2 in vduc_gx and g1 in mdo_gx:
        bezier(ax_rib, vduc_gx[g2], Y_VDUC + TH/2, mdo_gx[g1], Y_MDO - TH/2, '#A5B5C0', 0.13, 0.25, 2)
for g1, g2 in random.sample(vdar_vduc, min(2500, len(vdar_vduc))):
    if g1 in vdar_gx and g2 in vduc_gx:
        bezier(ax_rib, vdar_gx[g1], Y_VDAR + TH/2, vduc_gx[g2], Y_VDUC - TH/2, '#A5B5C0', 0.13, 0.25, 2)
    elif g2 in vdar_gx and g1 in vduc_gx:
        bezier(ax_rib, vdar_gx[g2], Y_VDAR + TH/2, vduc_gx[g1], Y_VDUC - TH/2, '#A5B5C0', 0.13, 0.25, 2)

n_almt_vm = n_almt_dv = 0
for g1, g2 in vduc_mdo:
    if g1 not in almt_ids or g2 not in almt_ids: continue
    if g1 in vduc_gx and g2 in mdo_gx: x1, x2 = vduc_gx[g1], mdo_gx[g2]
    elif g2 in vduc_gx and g1 in mdo_gx: x1, x2 = vduc_gx[g2], mdo_gx[g1]
    else: continue
    bezier(ax_rib, x1, Y_VDUC + TH/2, x2, Y_MDO - TH/2, '#C0392B', 0.85, 1.8, 8)
    n_almt_vm += 1
for g1, g2 in vdar_vduc:
    if g1 not in almt_ids or g2 not in almt_ids: continue
    if g1 in vdar_gx and g2 in vduc_gx: x1, x2 = vdar_gx[g1], vduc_gx[g2]
    elif g2 in vdar_gx and g1 in vduc_gx: x1, x2 = vdar_gx[g2], vduc_gx[g1]
    else: continue
    bezier(ax_rib, x1, Y_VDAR + TH/2, x2, Y_VDUC - TH/2, '#7B1FA2', 0.85, 1.8, 8)
    n_almt_dv += 1

def mark_almt(ax, gx, y, h):
    for gid in gx:
        if gid in almt_ids:
            x = gx[gid]
            ax.plot([x, x], [y - h/2 - 0.18, y - h/2 - 0.02], color='black', linewidth=1.3, zorder=12)
            ax.plot([x, x], [y + h/2 + 0.02, y + h/2 + 0.18], color='black', linewidth=1.3, zorder=12)

mark_almt(ax_rib, mdo_gx,  Y_MDO,  TH)
mark_almt(ax_rib, vduc_gx, Y_VDUC, TH)
mark_almt(ax_rib, vdar_gx, Y_VDAR, TH)

ax_rib.text(-0.05, 1.02, 'A', transform=ax_rib.transAxes, fontsize=22, fontweight='bold', va='top', ha='left')
legend_elements = [
    plt.Line2D([0], [0], color='#C0392B', linewidth=2.5, label=f'Vduc-Mdo ALMT syntenic pair (n={n_almt_vm})'),
    plt.Line2D([0], [0], color='#7B1FA2', linewidth=2.5, label=f'Vdar-Vduc ALMT syntenic pair (n={n_almt_dv})'),
    plt.Line2D([0], [0], color='#A5B5C0', linewidth=2,    label='Background syntenic blocks (subsampled n=2,500 per panel)'),
    plt.Line2D([0], [0], color='black', marker='|', markersize=12, linestyle='', label='ALMT gene position'),
]
ax_rib.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.07),
              ncol=2, fontsize=11, framealpha=0.95, handlelength=2.5)
ax_rib.axis('off')

mode_colors = {"dispersed":"#5DADE2","proximal":"#F39C12","tandem":"#9B59B6","WGD/segmental":"#E74C3C"}
sp_data = {
    "V. darrowii\n(n=8)":    {"dispersed":7, "WGD/segmental":1},
    "V. duclouxii\n(n=9)":   {"dispersed":6, "WGD/segmental":3},
    "M. domestica\n(n=26)":  {"proximal":2, "tandem":4, "WGD/segmental":20},
}
sp_baseline = {"V. darrowii\n(n=8)":17.9, "V. duclouxii\n(n=9)":13.8, "M. domestica\n(n=26)":37.8}
modes_order = ["dispersed","proximal","tandem","WGD/segmental"]
sp_pct = {sp: {m: d.get(m,0)/sum(d.values())*100 for m in modes_order} for sp, d in sp_data.items()}
species_list = list(sp_pct.keys())
x_b = np.arange(len(species_list)); bar_w = 0.45
bottom = np.zeros(len(species_list))
for mode in modes_order:
    vals = [sp_pct[sp][mode] for sp in species_list]
    if sum(vals) == 0: continue
    ax_dup.bar(x_b, vals, width=bar_w, bottom=bottom, color=mode_colors[mode],
               label=mode, edgecolor='white', linewidth=2.5)
    for i, v in enumerate(vals):
        if v > 4:
            raw_n = sp_data[species_list[i]].get(mode, 0)
            ax_dup.text(x_b[i], bottom[i] + v/2, f"{raw_n}",
                        ha='center', va='center', fontsize=12, fontweight='bold', color='white')
    bottom += vals
for i, sp in enumerate(species_list):
    by = sp_baseline[sp]
    ax_dup.plot([x_b[i] - bar_w/2 - 0.06, x_b[i] + bar_w/2 + 0.06], [by, by],
                color='#222222', linestyle=':', linewidth=2.5, zorder=10)
    # Compact label: just the percentage, placed on the line
    ax_dup.text(x_b[i] + bar_w/2 + 0.08, by, f"{by:.1f}%",
                fontsize=10, va='center', ha='left', color='#222222',
                fontweight='bold', zorder=11)
ax_dup.set_xticks(x_b); ax_dup.set_xticklabels(species_list, fontsize=12)
ax_dup.set_ylabel("% of ALMT genes", fontsize=12); ax_dup.set_ylim(0, 108)
ax_dup.text(-0.075, 1.06, 'B', transform=ax_dup.transAxes, fontsize=22, fontweight='bold', va='top', ha='left')
# Build duplication-mode legend with explicit baseline indicator
from matplotlib.lines import Line2D
mode_handles = [plt.Rectangle((0,0),1,1, facecolor=mode_colors[m],
                              edgecolor='white', linewidth=1.5) for m in modes_order]
mode_labels = list(modes_order)
# Add baseline line as the 5th entry
mode_handles.append(Line2D([0], [0], color='#222222', linestyle=':', linewidth=2.5))
mode_labels.append("Genome-wide WGD baseline (%)")
ax_dup.legend(mode_handles, mode_labels,
              loc='upper right', fontsize=10, framealpha=0.95,
              title="Duplication mode", title_fontsize=11,
              bbox_to_anchor=(1.0, 1.0))
ax_dup.grid(True, axis='y', alpha=0.3, linestyle=':')
ax_dup.spines['top'].set_visible(False); ax_dup.spines['right'].set_visible(False)
ax_dup.set_xlim(-0.5, 3.6)

plt.tight_layout()
for fmt in ("svg", "pdf", "png"):
    out = f"{fig_dir}/Figure4_synteny_duplication.{fmt}"
    plt.savefig(out, format=fmt, dpi=300 if fmt == "pdf" else 200, bbox_inches='tight')
    print(f"  wrote {out}")
plt.close()

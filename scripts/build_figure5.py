#!/usr/bin/env python3
"""Figure 5: PAML BEB plot. v5 + SVG + label fix.
Fix vs original: BEB-significant site labels now use Vduc1 REFERENCE position
(F249, E410, L411) consistent with manuscript text. Alignment column shown as
secondary annotation (col 325/487/488)."""
import os, re
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO

proj = os.path.expanduser("~/projects/vaccinium_almt")
paml = f"{proj}/results/p5_evolution/paml_og0008930"
fig_dir = f"{proj}/manuscript/figures"
os.makedirs(fig_dir, exist_ok=True)

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'E':-3.5,'Q':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2,'X':0,'*':0}

aln = {r.id: str(r.seq) for r in SeqIO.parse(f"{paml}/OG0008930.pep.aln", "fasta")}
ref_aligned = aln["Vduc1"]
ref_seq = ref_aligned.replace("-", "")
aln_len = len(ref_aligned)

W = 19; half = W // 2
kd_ref = []
for i in range(len(ref_seq)):
    win = ref_seq[max(0,i-half):min(len(ref_seq),i+half+1)]
    kd_ref.append(sum(KD.get(c,0) for c in win) / len(win))

THR, MIN_LEN = 0.8, 12
in_tm = False; start = None
tm_segs_ref = []
for i, v in enumerate(kd_ref):
    if v > THR and not in_tm: start = i; in_tm = True
    elif v <= THR and in_tm:
        if i - start >= MIN_LEN: tm_segs_ref.append((start, i))
        in_tm = False; start = None
if in_tm and len(kd_ref) - start >= MIN_LEN: tm_segs_ref.append((start, len(kd_ref)))

ref_to_aln = {}
rp = 0
for ap, aa in enumerate(ref_aligned):
    if aa != '-':
        ref_to_aln[rp] = ap; rp += 1
# Inverse: alignment column (0-indexed) -> Vduc1 1-indexed residue
aln_to_ref = {ap: rp + 1 for rp, ap in ref_to_aln.items()}

tm_segs_aln = []
for s, e in tm_segs_ref:
    if s in ref_to_aln and (e-1) in ref_to_aln:
        tm_segs_aln.append((ref_to_aln[s], ref_to_aln[e-1] + 1))

beb_sites = []
in_beb = False
with open(f"{paml}/modelA.out") as f:
    for line in f:
        if "Bayes Empirical Bayes" in line: in_beb = True; continue
        if in_beb:
            if "The grid" in line: break
            m = re.match(r'\s*(\d+)\s+(\S)\s+(\d+\.\d+)', line)
            if m: beb_sites.append((int(m.group(1)), m.group(2), float(m.group(3))))

fig = plt.figure(figsize=(15, 8.5))
ax_stats = fig.add_axes([0.08, 0.93, 0.88, 0.04])
ax_top   = fig.add_axes([0.08, 0.81, 0.88, 0.085])
ax_label = fig.add_axes([0.08, 0.73, 0.88, 0.05])
ax       = fig.add_axes([0.08, 0.09, 0.88, 0.59])

ax_stats.text(0.5, 0.5,
    "2$\\Delta$lnL = 40.65,  df = 1,  $\\it{p}$ = 9.1 x 10$^{-11}$    .    "
    "foreground $\\omega_2$ = 34.95 in 5.6% of sites    .    "
    "3 sites with BEB > 0.95",
    ha='center', va='center', fontsize=13)
ax_stats.axis('off')

ax_top.set_xlim(0, aln_len); ax_top.set_ylim(0, 1)
ax_top.axhline(y=0.5, xmin=0, xmax=1, color='#34495E', linewidth=2.8, zorder=1)
for i, (s, e) in enumerate(tm_segs_aln, 1):
    rect = Rectangle((s, 0.13), e-s, 0.74, facecolor='#E67E22', edgecolor='#7D3C00', linewidth=1.4, zorder=5)
    ax_top.add_patch(rect)
    ax_top.text((s+e)/2, 1.1, f'TM{i}', ha='center', va='bottom', fontsize=11, fontweight='bold', color='#7D3C00')
ax_top.text(-aln_len*0.013, 0.5, 'N', ha='right', va='center', fontsize=13, fontweight='bold')
ax_top.text(aln_len*1.013, 0.5, 'C', ha='left', va='center', fontsize=13, fontweight='bold')
ax_top.axis('off')

ax_label.set_xlim(0, aln_len); ax_label.set_ylim(0, 1)
n_term_end = tm_segs_aln[4][1] + 5 if len(tm_segs_aln) >= 5 else 280
ax_label.plot([0, n_term_end], [0.85, 0.85], color='#7D3C00', linewidth=1.8)
ax_label.plot([0, 0], [0.65, 0.85], color='#7D3C00', linewidth=1.8)
ax_label.plot([n_term_end, n_term_end], [0.65, 0.85], color='#7D3C00', linewidth=1.8)
ax_label.text(n_term_end/2, 0.4, 'N-terminal TM bundle (pore-forming)',
              ha='center', va='top', fontsize=12, color='#7D3C00', style='italic', fontweight='bold')
ax_label.plot([n_term_end+5, aln_len], [0.85, 0.85], color='#1F618D', linewidth=1.8)
ax_label.plot([n_term_end+5, n_term_end+5], [0.65, 0.85], color='#1F618D', linewidth=1.8)
ax_label.plot([aln_len, aln_len], [0.65, 0.85], color='#1F618D', linewidth=1.8)
ax_label.text((n_term_end+5+aln_len)/2, 0.4,
              'C-terminal cytoplasmic tail (regulatory)',
              ha='center', va='top', fontsize=12, color='#1F618D', style='italic', fontweight='bold')
ax_label.axis('off')

ax.set_xlim(0, aln_len); ax.set_ylim(0, 1.1)
for s, e in tm_segs_aln:
    ax.axvspan(s, e, color='#FEF5E7', alpha=0.6, zorder=0)
ax.axhline(y=0.95, color='#C0392B', linestyle='--', linewidth=1.2, alpha=0.75, zorder=2)
ax.axhline(y=0.50, color='gray', linestyle=':', linewidth=0.9, alpha=0.6, zorder=2)
ax.text(aln_len*0.998, 0.97, 'BEB > 0.95 (significant)', fontsize=11, color='#C0392B',
        ha='right', va='bottom', fontweight='bold')
ax.text(aln_len*0.998, 0.52, 'BEB > 0.5 (suggestive)', fontsize=11, color='gray',
        ha='right', va='bottom')

for pos, res, prob in beb_sites:
    if prob > 0.95: color, lw = '#C0392B', 2.2
    elif prob > 0.7: color, lw = '#E67E22', 1.3
    else: color, lw = '#888888', 0.8
    ax.plot([pos, pos], [0, prob], color=color, linewidth=lw, zorder=5, solid_capstyle='round')

sig_sites = sorted([s for s in beb_sites if s[2] > 0.95], key=lambda s: s[0])
label_positions = {325: (325, 1.06, 'center'), 487: (430, 1.06, 'right'), 488: (545, 1.06, 'left')}

print("\nBEB site mapping:")
for pos, res, prob in sig_sites:
    vduc1_pos = aln_to_ref.get(pos - 1, pos)   # PAML is 1-indexed, dict is 0-indexed
    print(f"  alignment col {pos} ({res}, p={prob:.3f}) -> Vduc1 {res}{vduc1_pos}")
    ax.scatter(pos, prob, s=110, c='#C0392B', edgecolor='black', linewidth=1.0, marker='o', zorder=10)
    lx, ly, ha_align = label_positions.get(pos, (pos, 1.06, 'center'))
    # 主标签用 Vduc1 reference 编号 (跟 manuscript 一致); col XXX 作次要
    ax.annotate(f"{res}{vduc1_pos}\n(col {pos})\n$\\it{{p}}$ = {prob:.3f}",
                xy=(pos, prob), xytext=(lx, ly),
                ha=ha_align, va='bottom',
                fontsize=11, fontweight='bold', color='#7B241C',
                arrowprops=dict(arrowstyle='-', color='#7B241C', linewidth=0.8, shrinkA=0, shrinkB=4))

ax.set_xlabel("Alignment column position (amino acid)", fontsize=13)
ax.set_ylabel("BEB posterior probability ($\\omega$ > 1)", fontsize=13)
ax.grid(True, axis='y', alpha=0.25, linestyle=':')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=11)

table_text = (
    "Site class       prop    bg w    fg w\n"
    "-------------------------------------\n"
    "  0 (conserved)  0.673   0.19    0.19\n"
    "  1 (neutral)    0.271   1.00    1.00\n"
    "  2a (selection) 0.040   0.19   34.95\n"
    "  2b (selection) 0.016   1.00   34.95"
)
ax.text(0.015, 0.97, table_text, transform=ax.transAxes,
        fontsize=10.5, family='monospace', verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.6', facecolor='#FDFEFE',
                  edgecolor='#34495E', linewidth=1.0, alpha=0.95))

for fmt in ("svg", "pdf", "png"):
    out = f"{fig_dir}/Figure5_paml_branchsite.{fmt}"
    plt.savefig(out, format=fmt, dpi=300 if fmt == "pdf" else 200, bbox_inches='tight')
    print(f"  wrote {out}")
plt.close()

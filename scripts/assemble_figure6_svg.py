#!/usr/bin/env python
"""
Figure 6 SVG: hybrid SVG with embedded PNG + vector labels/arrows/legend.
All text and annotations are editable in Inkscape/Illustrator.
"""
import os, math, base64
import numpy as np
from PIL import Image
from scipy.ndimage import label

P7 = "/home/us001/projects/vaccinium_almt/results/p7_structure"
OUTPUT_SIZE = 800
GAP = 30
MARGIN = 40

def find_red_clusters_raw(img_path, min_size_raw=500):
    img = np.array(Image.open(img_path).convert("RGB"))
    R, G, B = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
    is_red = (R > 120) & (G < 100) & (B < 100) & ((R - G) > 50)
    labeled, n = label(is_red)
    clusters = []
    for i in range(1, n+1):
        ys, xs = np.where(labeled == i)
        if len(ys) >= min_size_raw:
            clusters.append((int(xs.mean()), int(ys.mean()), len(ys)))
    return sorted(clusters, key=lambda c: -c[2])

def scale_to_800(cx, cy):
    return int(cx * 800 / 1600), int(cy * 800 / 1600)

def png_to_base64(png_path):
    """Encode PNG as base64 for SVG embedding."""
    with open(png_path, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")

# === Detect BEB site positions from raw PyMOL renders ===
clusters_B = find_red_clusters_raw(f"{P7}/panel_B_raw.png", min_size_raw=500)
clusters_C = find_red_clusters_raw(f"{P7}/panel_C_raw.png", min_size_raw=500)
clusters_D = find_red_clusters_raw(f"{P7}/panel_D_raw.png", min_size_raw=500)

# B: identify isolated F249 vs E410/L411 pair
top3_B = clusters_B[:3]
def min_dist(c, others):
    return min(((c[0]-o[0])**2 + (c[1]-o[1])**2)**0.5 for o in others if o != c)
isolation = sorted([(c, min_dist(c, top3_B)) for c in top3_B], key=lambda x: -x[1])
f249_cluster = isolation[0][0]
el_clusters = [c for c in top3_B if c != f249_cluster]
f249_xy = scale_to_800(f249_cluster[0], f249_cluster[1])
el_avg = (int(np.mean([c[0] for c in el_clusters])), int(np.mean([c[1] for c in el_clusters])))
el_xy = scale_to_800(*el_avg)

# C: single cluster (the BEB sites group viewed from rotated angle)
c_xy = scale_to_800(clusters_C[0][0], clusters_C[0][1]) if clusters_C else (400, 400)

# D: two clusters (E410 + L411)
big_D = [c for c in clusters_D if c[2] > 1000][:2]
big_D_sortx = sorted(big_D, key=lambda c: c[0])
d_left = scale_to_800(big_D_sortx[0][0], big_D_sortx[0][1])
d_right = scale_to_800(big_D_sortx[1][0], big_D_sortx[1][1])

# === Canvas dimensions ===
W = OUTPUT_SIZE*2 + GAP + MARGIN*2
H = OUTPUT_SIZE*2 + GAP + MARGIN*2

# Panel positions (top-left corner of each panel in canvas)
A_pos = (MARGIN, MARGIN)
B_pos = (MARGIN + OUTPUT_SIZE + GAP, MARGIN)
C_pos = (MARGIN, MARGIN + OUTPUT_SIZE + GAP)
D_pos = (MARGIN + OUTPUT_SIZE + GAP, MARGIN + OUTPUT_SIZE + GAP)

# Encode panel PNGs
print("Encoding PNG panels as base64...")
A_b64 = png_to_base64(f"{P7}/panel_A_raw.png")
B_b64 = png_to_base64(f"{P7}/panel_B_raw.png")
C_b64 = png_to_base64(f"{P7}/panel_C_raw.png")
D_b64 = png_to_base64(f"{P7}/panel_D_raw.png")

# === Build SVG ===
# AlphaFold colors
AF_VHIGH = "rgb(13,87,211)"
AF_HIGH  = "rgb(106,203,240)"
AF_LOW   = "rgb(254,217,54)"
AF_VLOW  = "rgb(253,125,77)"

# Helper: arrow path
def arrow_marker():
    return '''<defs>
    <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto" markerUnits="strokeWidth">
        <polygon points="0 0, 10 3.5, 0 7" fill="#282828"/>
    </marker>
</defs>'''

# Helper: text with white outline (using paint-order)
def label_text(x, y, text, size=24, anchor="start"):
    return f'''<text x="{x}" y="{y}" font-family="DejaVu Sans, Arial, sans-serif" 
        font-size="{size}" font-weight="bold" 
        fill="black" stroke="white" stroke-width="3" 
        paint-order="stroke fill" text-anchor="{anchor}">{text}</text>'''

# Helper: arrow line
def arrow(x1, y1, x2, y2):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#282828" stroke-width="3" marker-end="url(#arrowhead)"/>'

# Embed each panel image
def embed_panel(x, y, b64_data, panel_letter):
    return f'''<image x="{x}" y="{y}" width="{OUTPUT_SIZE}" height="{OUTPUT_SIZE}" 
        xlink:href="data:image/png;base64,{b64_data}" preserveAspectRatio="xMidYMid meet"/>
    <text x="{x+25}" y="{y+70}" font-family="DejaVu Sans, Arial, sans-serif" 
        font-size="60" font-weight="bold" fill="black">{panel_letter}</text>'''

# Compute label/arrow positions in CANVAS coordinates (panel offset + position-within-panel)
def to_canvas(panel_pos, xy):
    return (panel_pos[0] + xy[0], panel_pos[1] + xy[1])

# Panel B label positions (within-panel coords)
b_f249_canvas = to_canvas(B_pos, f249_xy)
b_el_canvas = to_canvas(B_pos, el_xy)

b_f249_label_pos = (max(B_pos[0]+20, b_f249_canvas[0] - 180), max(B_pos[1]+20, b_f249_canvas[1] - 30))
b_f249_arrow_start = (b_f249_label_pos[0] + 95, b_f249_label_pos[1] + 12)
b_f249_arrow_end = (b_f249_canvas[0] - 30, b_f249_canvas[1])

b_el_label_pos = (b_el_canvas[0] + 60, b_el_canvas[1] + 30)
b_el_arrow_start = (b_el_label_pos[0] - 8, b_el_label_pos[1] - 8)
b_el_arrow_end = (b_el_canvas[0] + 20, b_el_canvas[1])

# Panel C label
c_canvas = to_canvas(C_pos, c_xy)
c_label_pos = (c_canvas[0] + 80, c_canvas[1] - 30)
c_arrow_start = (c_label_pos[0] - 8, c_label_pos[1] + 8)
c_arrow_end = (c_canvas[0] + 20, c_canvas[1])

# Panel D labels  
d_left_canvas = to_canvas(D_pos, d_left)
d_right_canvas = to_canvas(D_pos, d_right)
d_e410_label_pos = (d_left_canvas[0] - 150, d_left_canvas[1] - 80)
d_e410_arrow_start = (d_e410_label_pos[0] + 30, d_e410_label_pos[1] + 12)
d_e410_arrow_end = (d_left_canvas[0] - 10, d_left_canvas[1] - 10)

d_l411_label_pos = (d_right_canvas[0] + 80, d_right_canvas[1] + 80)
d_l411_arrow_start = (d_l411_label_pos[0], d_l411_label_pos[1] - 8)
d_l411_arrow_end = (d_right_canvas[0] + 10, d_right_canvas[1] + 5)

# Panel A legend
legend_x = A_pos[0] + 30
legend_y = A_pos[1] + 600
legend_items = [
    (AF_VHIGH, "Very high (≥90)"),
    (AF_HIGH,  "Confident (70-90)"),
    (AF_LOW,   "Low (50-70)"),
    (AF_VLOW,  "Very low (<50)"),
]
legend_svg = f'''<g id="panel_A_legend">
    <rect x="{legend_x-5}" y="{legend_y-5}" width="230" height="175" 
          fill="white" stroke="rgb(180,180,180)" stroke-width="2"/>
    <text x="{legend_x+5}" y="{legend_y+22}" font-family="DejaVu Sans, Arial, sans-serif" 
          font-size="20" font-weight="bold" fill="black">Model confidence</text>'''
for i, (col, txt) in enumerate(legend_items):
    yy = legend_y + 40 + i*33
    legend_svg += f'''
    <rect x="{legend_x+10}" y="{yy}" width="30" height="28" fill="{col}" stroke="black" stroke-width="1"/>
    <text x="{legend_x+50}" y="{yy+20}" font-family="DejaVu Sans, Arial, sans-serif" 
          font-size="18" fill="black">{txt}</text>'''
legend_svg += '\n</g>'

# === Assemble SVG ===
svg = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg" 
     xmlns:xlink="http://www.w3.org/1999/xlink"
     width="{W}" height="{H}" viewBox="0 0 {W} {H}">

<title>Figure 6 - Vduc1 ALMT BEB Sites</title>

{arrow_marker()}

<rect width="{W}" height="{H}" fill="white"/>

<!-- Panel A: pLDDT -->
<g id="panel_A">
{embed_panel(A_pos[0], A_pos[1], A_b64, "A")}
{legend_svg}
</g>

<!-- Panel B: BEB sites overview -->
<g id="panel_B">
{embed_panel(B_pos[0], B_pos[1], B_b64, "B")}

<!-- F249 label + arrow -->
{arrow(b_f249_arrow_start[0], b_f249_arrow_start[1], b_f249_arrow_end[0], b_f249_arrow_end[1])}
{label_text(b_f249_label_pos[0], b_f249_label_pos[1] + 28, "F249", size=32)}

<!-- E410/L411 label + arrow -->
{arrow(b_el_arrow_start[0], b_el_arrow_start[1], b_el_arrow_end[0], b_el_arrow_end[1])}
{label_text(b_el_label_pos[0], b_el_label_pos[1] + 28, "E410/L411", size=32)}
</g>

<!-- Panel C: rotated view -->
<g id="panel_C">
{embed_panel(C_pos[0], C_pos[1], C_b64, "C")}
{arrow(c_arrow_start[0], c_arrow_start[1], c_arrow_end[0], c_arrow_end[1])}
{label_text(c_label_pos[0], c_label_pos[1] + 28, "BEB sites", size=32)}
</g>

<!-- Panel D: C-terminal close-up -->
<g id="panel_D">
{embed_panel(D_pos[0], D_pos[1], D_b64, "D")}

<!-- E410 -->
{arrow(d_e410_arrow_start[0], d_e410_arrow_start[1], d_e410_arrow_end[0], d_e410_arrow_end[1])}
{label_text(d_e410_label_pos[0], d_e410_label_pos[1] + 28, "E410", size=32)}

<!-- L411 -->
{arrow(d_l411_arrow_start[0], d_l411_arrow_start[1], d_l411_arrow_end[0], d_l411_arrow_end[1])}
{label_text(d_l411_label_pos[0], d_l411_label_pos[1] + 28, "L411", size=32)}
</g>

</svg>
'''

out = f"{P7}/Figure6_final.svg"
with open(out, "w") as f:
    f.write(svg)

print(f"\n✓ SVG saved: {out}")
print(f"  Size: {os.path.getsize(out)/1024:.0f} KB")
print(f"  Canvas: {W} x {H}")
print(f"\nDetected positions (canvas coords):")
print(f"  Panel B F249 sphere: {b_f249_canvas}")
print(f"  Panel B E410/L411 sphere: {b_el_canvas}")
print(f"  Panel C BEB cluster: {c_canvas}")
print(f"  Panel D E410: {d_left_canvas}")
print(f"  Panel D L411: {d_right_canvas}")

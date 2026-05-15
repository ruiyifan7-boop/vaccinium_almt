# Figure 6 v6 - Fix A/B occlusion: side-view orient + transparent backdrop for B
load Vduc1_OG0008930_ecbc1_unrelaxed_alphafold2_ptm_model_1_seed_000.pdb, vduc1

# === GLOBAL SETUP ===
bg_color white
set opaque_background, on
set ray_opaque_background, on
set ray_shadows, 0
set ray_trace_mode, 1
set ray_trace_color, gray30
set antialias, 2
set ambient, 0.3
set specular, 0.2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# === Custom colours ===
set_color af_vhigh, [0.051, 0.341, 0.827]
set_color af_high,  [0.416, 0.796, 0.945]
set_color af_low,   [0.996, 0.851, 0.212]
set_color af_vlow,  [0.992, 0.490, 0.302]
set_color tm_color, [0.42, 0.55, 0.80]
set_color loop_color, [0.85, 0.78, 0.55]
set_color ctail_color, [0.95, 0.65, 0.55]
set_color beb_red, [0.85, 0.15, 0.15]

# ====================================================================
# PANEL A: pLDDT discrete coloring + SIDE VIEW (TM helices horizontal)
# ====================================================================
hide everything
show cartoon

select vlow_residues, byres (name CA and b<50)
select low_or_below, byres (name CA and b<70)
select low_residues, low_or_below and not vlow_residues
select high_or_below, byres (name CA and b<90)
select high_residues, high_or_below and not low_or_below
select vhigh_residues, byres (name CA and not high_or_below)
color af_vlow, vlow_residues
color af_low, low_residues
color af_high, high_residues
color af_vhigh, vhigh_residues
set cartoon_loop_radius, 0.15, vlow_residues

# Side view: orient + rotate so TM bundle lies horizontally
# orient = align longest principal axis with x-axis (default)
# Without extra turns, this gives a side-view of the membrane plane.
orient
# Slight tilt to add depth (avoid pure flat side view)
turn y, 15
turn x, 5
zoom all, 3

ray 1600, 1600
png figure6_panelA_pLDDT.png, dpi=300
print ">>> Panel A saved (side view, no occlusion)"

# ====================================================================
# PANEL B: BEB sites — transparent backdrop, large spheres
# ====================================================================
hide everything
show cartoon
color gray70, all

select TM_all, resi 49-66 or resi 75-93 or resi 103-143 or resi 186-204 or resi 383-402
color tm_color, TM_all
select TM4_TM5_loop, resi 205-382
color loop_color, TM4_TM5_loop
select C_tail, resi 403-420
color ctail_color, C_tail
set cartoon_loop_radius, 0.15, vlow_residues

# KEY FIX: make entire cartoon semi-transparent so BEB spheres show through
set cartoon_transparency, 0.55, all

# BEB sites (opaque, large, black-edged)
select beb_F249, resi 249
select beb_E410, resi 410
select beb_L411, resi 411
select beb_all, beb_F249 or beb_E410 or beb_L411
show spheres, beb_all and name CA
color beb_red, beb_all
set sphere_scale, 3.0, beb_all
# Black outline for spheres (set sphere outline)
set ray_trace_disabled_shadows, 1
set sphere_quality, 3

# View: orient on BEB sites, then back off to show full protein
orient beb_all
zoom vduc1, 5
turn y, 20

python
f249_xyz = cmd.get_atom_coords("resi 249 and name CA")
e410_xyz = cmd.get_atom_coords("resi 410 and name CA")
l411_xyz = cmd.get_atom_coords("resi 411 and name CA")
cmd.pseudoatom("lbl_F249", pos=[f249_xyz[0]-25, f249_xyz[1]+25, f249_xyz[2]])
cmd.pseudoatom("lbl_E410", pos=[e410_xyz[0]+5, e410_xyz[1]-30, e410_xyz[2]])
cmd.pseudoatom("lbl_L411", pos=[l411_xyz[0]+30, l411_xyz[1]+5, l411_xyz[2]])
python end

hide everything, lbl_*
label lbl_F249, "F249"
label lbl_E410, "E410"
label lbl_L411, "L411"
set label_size, 32
set label_color, black
set label_font_id, 7
set label_outline_color, white

ray 1800, 1600
png figure6_panelB_BEBsites.png, dpi=300
print ">>> Panel B saved (transparent backdrop)"

# Reset transparency for subsequent panels
set cartoon_transparency, 0, all

# ====================================================================
# PANEL C: C-tail close-up (unchanged from v5)
# ====================================================================
hide everything
hide labels
delete lbl_*
select c_region, resi 395-420
show cartoon, c_region
color ctail_color, c_region
select kevdel, resi 405-411
color salmon, kevdel
show sticks, kevdel and not (name N+C+O+H*)
color gray40, kevdel and elem C
select beb_pair, resi 410 or resi 411
show sticks, beb_pair and not (name N+C+O+H*)
color beb_red, beb_pair and elem C
show spheres, beb_pair and name CA
set sphere_scale, 0.8, beb_pair
orient c_region

python
e410_xyz = cmd.get_atom_coords("resi 410 and name CA")
l411_xyz = cmd.get_atom_coords("resi 411 and name CA")
cmd.pseudoatom("cc_E410", pos=[e410_xyz[0]-5, e410_xyz[1]-12, e410_xyz[2]])
cmd.pseudoatom("cc_L411", pos=[l411_xyz[0]+5, l411_xyz[1]+12, l411_xyz[2]])
python end

hide everything, cc_*
label cc_E410, "E410"
label cc_L411, "L411"
set label_size, 40
set label_color, black
set label_font_id, 7
set label_outline_color, white
zoom c_region or cc_*, 3
ray 1600, 1200
png figure6_panelC_Ctail.png, dpi=300
print ">>> Panel C saved"

# ====================================================================
# PANEL D: F249 close-up (unchanged)
# ====================================================================
hide everything
hide labels
delete cc_*
select f249_region, resi 230-270
show cartoon, f249_region
color loop_color, f249_region
select f249, resi 249
show spheres, f249 and name CA
color beb_red, f249
set sphere_scale, 1.0, f249
show sticks, f249 and not (name N+C+O+H*)
color beb_red, f249 and elem C
orient f249_region

python
f249_xyz = cmd.get_atom_coords("resi 249 and name CA")
cmd.pseudoatom("dd_F249", pos=[f249_xyz[0], f249_xyz[1]+15, f249_xyz[2]])
python end

hide everything, dd_F249
label dd_F249, "F249"
set label_size, 40
set label_color, black
set label_font_id, 7
set label_outline_color, white
zoom f249_region or dd_F249, 3
ray 1600, 1200
png figure6_panelD_F249.png, dpi=300
print ">>> Panel D saved"
delete dd_F249

print ""
print "=== Figure 6 v6 done ==="

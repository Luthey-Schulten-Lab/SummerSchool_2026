# =============================================================================
#  representations.tcl  --  4DWCM (STC-QCB summer school 2026)
#
#  Visual styling for the two molecules loaded by load_and_sync.tcl.
#  Expects these globals (set by load_and_sync.tcl):
#      dna_mol   -- molid of the LAMMPS chromosome
#      lm_mol    -- molid of the LM lattice  (optional)
#
#  Standalone use:  set the molids by hand, then source this file, e.g.
#      set dna_mol [molinfo top] ; source representations.tcl
# =============================================================================

# ----------------------------- global display -------------------------------
color Display Background white
color scale method viridis
color scale min 0.0
color scale midpoint 0.5
color scale max 1.0
display projection Orthographic
display height     4.0          ;# Display menu: Screen Hgt
display distance   -2.0         ;# Display menu: Screen Dist

# Depth-fog look carried over from movie.vmd (save_state).
display depthcue         on
display cuemode          Linear
display cuestart         0.5
display cueend           10.0
display cuedensity       0.32
display shadows          off
display ambientocclusion off


# ------------------- DNA (LAMMPS) per-bead-type representations --------------
#  Bead types are carried in vy (= c_type_track, set via LAMMPSREMAPFIELDS):
#    1 bdry   2 ribo   3 mono(DNA)   4 ori   5 ter   6 fork   7 anchor   8 hinge
proc make_dna_representations {molid} {
    # --- bead radii (Angstrom) ---
    set DNA_rad    13.0
    set a_special  2.0
    set ori_rad    [expr {$a_special*$DNA_rad}]
    set ter_rad    [expr {$a_special*$DNA_rad}]
    set fork_rad   [expr {$a_special*$DNA_rad}]
    set hinge_rad  [expr {1.5*$DNA_rad}]
    set anchor_rad [expr {1.5*$DNA_rad}]
    set ribo_rad   70.0
    set bdry_rad   65.0

    # --- selections (text only; selupdate keeps them live as the cell grows) ---
    set s_DNA    "vy>2"
    set s_DNA2   "vy>2 and index < 54338"
    set s_ori    "vy==4"
    set s_ter    "vy==5"
    set s_fork   "vy==6"
    set s_anchor "vy==7"
    set s_hinge  "vy==8"
    set s_ribo   "vy==2"
    set s_bdry   "vy==1"

    # --- clear any existing reps on this molecule ---
    for {set r [expr {[molinfo $molid get numreps]-1}]} {$r >= 0} {incr r -1} {
        mol delrep $r $molid
    }

    # --- one VDW rep per type:  {selection radius resolution material colorID show} ---
    #  Colours and visibility carried over from movie.vmd (GUI-tuned):
    #    * DNA recoloured red (1), first arm recoloured blue (0).
    #    * fork/anchor/hinge/ribo/bdry hidden -- ribosomes + membrane are drawn
    #      from the LM lattice instead, so we don't double-draw them here.
    set rep 0
    foreach {sel radius res mat cid show} [list \
            $s_DNA    $DNA_rad    12.0 AOChalky     1 on  \
            $s_DNA2   13.1        12.0 AOChalky     0 on  \
            $s_ori    $ori_rad    20.0 AOChalky     1 on  \
            $s_ter    $ter_rad    20.0 AOChalky     3 on  \
            $s_fork   $fork_rad   20.0 AOChalky    27 off \
            $s_anchor $anchor_rad 12.0 AOChalky    16 off \
            $s_hinge  $hinge_rad  12.0 AOChalky     8 off \
            $s_ribo   $ribo_rad   20.0 AOChalky    13 off \
            $s_bdry   $bdry_rad   20.0 Transparent  2 off ] {
        mol addrep      $molid
        mol modstyle    $rep $molid VDW $radius $res
        mol modmaterial $rep $molid $mat
        mol modselect   $rep $molid $sel
        mol modcolor    $rep $molid ColorID $cid
        mol selupdate   $rep $molid on
        mol colupdate   $rep $molid on
        mol showrep     $molid $rep $show
        incr rep
    }
    puts "  DNA representations: created $rep reps on mol $molid"
}


# --------------------------- LM lattice representations ----------------------
#  Carried over from movie.vmd (GUI-tuned).  Three reps:
#    0  membrane as LatticeCubes, near half cut away (y<3200) to see inside
#    1  ribosome particles (species 'ribosomeP') as CPK
#    2  individually-named ribosomes RB_#### (exactly four digits) as CPK
#  selupdate is forced ON for all three: RDME particles move/react every frame
#  and the cell grows, so the selections must stay live (movie.vmd had it off
#  for reps 1-2 -- flip back to off if you want a frozen first-frame selection).
proc make_lm_representations {molid} {
    # --- clear any existing reps on this molecule ---
    for {set r [expr {[molinfo $molid get numreps]-1}]} {$r >= 0} {incr r -1} {
        mol delrep $r $molid
    }

    # rep 0: membrane lattice cubes, cutaway lower half
    mol addrep      $molid
    mol modstyle    0 $molid LatticeCubes 1.0
    mol modcolor    0 $molid ColorID 7
    mol modselect   0 $molid {name membrane and y < 3200}
    mol modmaterial 0 $molid Opaque
    mol selupdate   0 $molid on
    mol colupdate   0 $molid on

    # rep 1: ribosome species 'ribosomeP'
    mol addrep      $molid
    mol modstyle    1 $molid CPK 20.0 0.3 12.0 12.0
    mol modcolor    1 $molid ColorID 11
    mol modselect   1 $molid {name ribosomeP}
    mol modmaterial 1 $molid Opaque
    mol selupdate   1 $molid on
    mol colupdate   1 $molid on

    # rep 2: individually-named ribosomes RB_#### (exactly four digits)
    mol addrep      $molid
    mol modstyle    2 $molid CPK 20.0 0.3 12.0 12.0
    mol modcolor    2 $molid ColorID 4
    mol modselect   2 $molid {name "RB_[0-9][0-9][0-9][0-9]"}
    mol modmaterial 2 $molid Opaque
    mol selupdate   2 $molid on
    mol colupdate   2 $molid on

    puts "  LM representations: created 3 reps on mol $molid"
}


# ---------------------------------- apply -----------------------------------
if {[info exists dna_mol]} {
    make_dna_representations $dna_mol
} else {
    puts "representations.tcl: \$dna_mol not set; skipping DNA reps."
}
if {[info exists lm_mol]} {
    make_lm_representations $lm_mol
}

# ---------------------------------- camera ----------------------------------
# Saved camera angle from movie.vmd: {center_matrix rotate_matrix scale_matrix
# global_matrix}.  Applied to both molecules so the scene keeps the tuned view
# (replaces `display resetview`, which would reset it to front-on).
#
# The rotate_matrix is the movie.vmd one FLIPPED 180 deg about the screen-X axis
# (Y and Z rows negated => Y-diagonal -0.907 -> +0.907), so model +Y now points
# UP.  This puts the physical lower half (membrane `y < 3200`) at the BOTTOM of
# the window, with the cut-away (open) side facing the camera.  To go back to the
# original movie.vmd orientation, re-negate the Y and Z rows below.
set _saved_view {{{1 0 0 -6492.51} {0 1 0 -3257.87} {0 0 1 -3251.39} {0 0 0 1}} {{0.999738 -0.0184609 -0.0139979 0} {0.0108502 0.906948 -0.421102 0} {0.0204689 0.420838 0.906911 0} {0 0 0 1}} {{0.000242125 0 0 0} {0 0.000242125 0 0} {0 0 0.000242125 0} {0 0 0 1}} {{1 0 0 0.00999999} {0 1 0 0.08} {0 0 1 0} {0 0 0 1}}}
foreach _m {lm_mol dna_mol} {
    if {[info exists $_m]} {
        molinfo [set $_m] set {center_matrix rotate_matrix scale_matrix global_matrix} $_saved_view
    }
}

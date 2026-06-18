# =============================================================================
#  load_and_sync.tcl   --   4DWCM (STC-QCB summer school 2026)
#
#  Loads and TIME-SYNCHRONIZES the two trajectories of one 4DWCM run so they
#  play frame-for-frame in lockstep, then (optionally) sources representations.tcl.
#
#  Replaces the manual recipe in note_andrew plus the three helper scripts:
#      translate_and_rotate.tcl + remove_bdryGrow_frames.tcl + keep_every_four.tcl
#
#  Background (whole simulation = 2 h = 7200 s of cell-cycle time):
#    * LM  (RDME lattice, MinCell.lm)      : 1 frame / second            (native)
#    * DNA (LAMMPS chromosome.lammpstrj)   : ~1 frame / 4 s AFTER we strip the
#      boundary-growth ("BdryGrow") frames; 4 s is its finest spacing.
#
#  Both are reduced to a COMMON cadence `lm_freq` (seconds, multiple of 4):
#      LM  : keep 1 frame every  lm_freq        seconds        (stride on a 1 s traj)
#      DNA : keep 1 frame every (lm_freq / 4)   cleaned frames (4 s/frame traj)
#  =>  frame k of LM  <->  frame k of DNA  <->  cell time  k * lm_freq  seconds.
#
#  Usage (VMD TkConsole, cwd = render_scripts):   source load_and_sync.tcl
# =============================================================================

# ------------------------------- PARAMETERS ----------------------------------
# Common render cadence, in seconds. Integer, multiple of 4, >= 4.
#   4 -> finest (LM and DNA both 1:1).  8, 12, ... -> coarser / lighter movies.
set lm_freq    60

# Total cell-cycle time to render, in seconds. Integer, multiple of 4, <= 7200,
# and should be a multiple of lm_freq.  (DNA cadence is a clean 4 s/frame in the
# middle of the cycle but drifts during the final period -- stop short of 7200 if
# you need perfect sync all the way to the end.)
set duration   7200

# Input files for this run -- edit per run.
set run_dir    ../trajectory/Mar31_1/
set lm_file    ../trajectory/Mar31_1/MinCell.lm
set dna_file   ../trajectory/Mar31_1/DNA/chromosome.lammpstrj

# Apply representations.tcl automatically after loading?  (1 = yes, 0 = no)
# Kept 0: load + sync only.  Apply visuals yourself with `source representations.tcl`.
set apply_reps 0
# -----------------------------------------------------------------------------


# ------------------------------- VALIDATION ----------------------------------
proc _check_mult4 {name val} {
    if {![string is integer -strict $val] || $val < 4 || $val % 4 != 0} {
        error "$name must be an integer multiple of 4 and >= 4 (got '$val')"
    }
}
_check_mult4 lm_freq  $lm_freq
_check_mult4 duration $duration
if {$duration % $lm_freq != 0} {
    puts "WARNING: duration ($duration s) is not a multiple of lm_freq ($lm_freq s); final partial step dropped."
}
set dna_stride [expr {$lm_freq / 4}]   ;# cleaned-DNA frames per rendered frame
puts "load_and_sync: cadence ${lm_freq}s/frame, duration ${duration}s -> [expr {$duration/$lm_freq + 1}] frames target (dna_stride $dna_stride)"


# ------------------------------- PROCEDURES ----------------------------------

# Remove LAMMPS boundary-growth frames: runs of frames that have NO anchor
# (vy==7) and NO hinge (vy==8) are collapsed to a single frame.  Faithful to
# remove_bdryGrow_frames.tcl (loop start fixed to the last valid frame index).
proc remove_bdrygrow_frames {molid} {
    set batch {}
    for {set i [expr {[molinfo $molid get numframes]-1}]} {$i >= 0} {incr i -1} {
        set anchor [atomselect $molid "vy==7" frame $i]
        set hinge  [atomselect $molid "vy==8" frame $i]
        if {[$anchor num]==0 && [$hinge num]==0} {
            lappend batch $i
        } elseif {[llength $batch] > 0} {
            set keep [lindex $batch end]
            foreach f $batch { if {$f != $keep} { animate delete beg $f end $f } }
            set batch {}
        }
        $anchor delete
        $hinge delete
    }
    if {[llength $batch] > 0} {
        set keep [lindex $batch end]
        foreach f $batch { if {$f != $keep} { animate delete beg $f end $f } }
    }
}

# Keep frames whose 0-based index is a multiple of `stride` and <= `last_idx`;
# delete the rest.  Afterwards frame j corresponds to original frame j*stride.
proc decimate_frames {molid stride last_idx} {
    for {set i [expr {[molinfo $molid get numframes]-1}]} {$i >= 0} {incr i -1} {
        if {$i > $last_idx || ($i % $stride) != 0} {
            animate delete beg $i end $i
        }
    }
}

# Move the DNA molecule into the LM lattice frame: swap X<->Z, then translate to
# the LM box centre.  (From translate_and_rotate.tcl lines 175-205.)
proc align_dna_to_lm {molid} {
    set all  [atomselect $molid all]
    set swap [list {0 0 1 0} {0 1 0 0} {1 0 0 0} {0 0 0 1}]
    set n    [molinfo $molid get numframes]
    for {set i 0} {$i < $n} {incr i} {
        $all frame $i
        $all move $swap
        $all moveby {6500.0 3250.0 3250.0}
    }
    $all delete
}


# --------------------------------- LOAD LM -----------------------------------
# LM is 1 frame/s, so frame index == seconds.  Load only [0 .. duration] at a
# stride of lm_freq so we never hold the whole 7201-frame trajectory in memory.
puts "Loading LM:  $lm_file   (last=$duration step=$lm_freq) ..."
set lm_mol [mol new $lm_file type LM first 0 last $duration step $lm_freq waitfor all]
mol rename $lm_mol "LM_lattice"
puts "  LM frames loaded: [molinfo $lm_mol get numframes]"


# --------------------------------- LOAD DNA ----------------------------------
# LAMMPS plugin env vars MUST be set before the load.
set env(LAMMPSDUMMYPOS)    "0.0,0.0,0.0"
set env(LAMMPSMAXATOMS)    200000
set env(LAMMPSREMAPFIELDS) {vx=c_id_track,vy=c_type_track}

# `step 1` is required (note_andrew) or only one frame loads.  ALL DNA frames
# must be loaded: BdryGrow frames are interleaved and cannot be strided out.
puts "Loading DNA: $dna_file   (all frames) ..."
set dna_mol [mol new $dna_file type lammpstrj first 0 last -1 step 1 waitfor all]
mol rename $dna_mol "DNA_chromosome"
puts "  DNA frames loaded: [molinfo $dna_mol get numframes]"


# ------------------------- CLEAN + SYNC THE DNA TRAJ -------------------------
puts "Cleaning DNA: removing boundary-growth frames ..."
remove_bdrygrow_frames $dna_mol
puts "  DNA frames after BdryGrow removal: [molinfo $dna_mol get numframes]  (expect ~1801 for the full run)"

# Reduce cleaned DNA (4 s/frame) to the common cadence and requested duration.
decimate_frames $dna_mol $dna_stride [expr {$duration/4}]
puts "  DNA frames after sync/decimation:  [molinfo $dna_mol get numframes]"

# Finally, line the DNA up in space with the LM lattice.
align_dna_to_lm $dna_mol


# ------------------------------- SANITY CHECK --------------------------------
set nlm  [molinfo $lm_mol  get numframes]
set ndna [molinfo $dna_mol get numframes]
puts "----------------------------------------------------------------------"
puts " SYNC SUMMARY:  LM frames = $nlm   DNA frames = $ndna   (cadence ${lm_freq}s)"
if {$nlm != $ndna} {
    puts " WARNING: frame counts differ by [expr {abs($nlm-$ndna)}] -> movie will drift."
    puts "          Lower `duration` to stay in the clean-4s region, or re-check"
    puts "          the DNA BdryGrow cleanup."
} else {
    puts " OK: trajectories are frame-synced  (frame k <-> t = k*${lm_freq}s)."
}
puts "----------------------------------------------------------------------"


# ------------------------------ REPRESENTATIONS ------------------------------
# Visuals live in representations.tcl so this file stays purely load + sync.
# That file uses the globals set above: $lm_mol and $dna_mol.
if {$apply_reps} {
    set _repfile [file join [file dirname [info script]] representations.tcl]
    if {[file exists $_repfile]} {
        puts "Applying representations from $_repfile ..."
        source $_repfile
    } else {
        puts "NOTE: $_repfile not found; skipping representations."
    }
}

display resetview

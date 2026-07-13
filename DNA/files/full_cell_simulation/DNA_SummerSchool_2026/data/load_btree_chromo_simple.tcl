## Simple VMD setup for modify_lammpstrj.py output (mother / left / right daughters).
## Run after: python3 modify_lammpstrj.py
##
## c_type_track in modified.lammpstrj:
##   3  = mother DNA (gray)
##  13  = left daughter (lime)
##  14  = right daughter / DNA_new (magenta)
##   7  = anchor / SMC1 (black)
##   8  = hinge / SMC2 (white)

set xd 0.0
set yd 0.0
set zd 0.0
set env(LAMMPSDUMMYPOS) {$xd,$yd,$zd}

set Nmax 200000
set env(LAMMPSMAXATOMS) $Nmax
set env(LAMMPSREMAPFIELDS) {vx=c_id_track,vy=c_type_track}

color Display {Background} white
display projection Orthographic

set DNA_rad 13.0
set a_special 2.0
set ori_rad [expr $a_special * $DNA_rad]
set ter_rad [expr $a_special * $DNA_rad]
set fork_rad [expr $a_special * $DNA_rad]
set hinge_rad [expr 1.5 * $DNA_rad]
set anchor_rad [expr 1.5 * $DNA_rad]
set ribo_rad 70.0
set bdry_rad 32.5

set molecule_file modified.lammpstrj
mol new $molecule_file waitfor all type lammpstrj step 1
set mol_id [molinfo top]

set bdry [atomselect $mol_id "vy==1" frame now]
set ribo [atomselect $mol_id "vy==2" frame now]
set DNA_mother [atomselect $mol_id "vy==3" frame now]
set DNA_left [atomselect $mol_id "vy==13" frame now]
set DNA_right [atomselect $mol_id "vy==14" frame now]
set ori [atomselect $mol_id "vy==4" frame now]
set ter [atomselect $mol_id "vy==5" frame now]
set fork [atomselect $mol_id "vy==6" frame now]
set anchor [atomselect $mol_id "vy==7" frame now]
set hinge [atomselect $mol_id "vy==8" frame now]

mol modstyle 0 $mol_id Points
mol delrep 0 $mol_id

set rep_count 0

proc add_vdw_rep {mol_id rep_count selection radius color} {
    mol addrep $mol_id
    set rep $rep_count
    mol modstyle $rep $mol_id VDW $radius 12.0
    mol modmaterial $rep $mol_id AOChalky
    mol modselect $rep $mol_id $selection
    mol selupdate $rep $mol_id on
    mol modcolor $rep $mol_id ColorID $color
    return [expr $rep_count + 1]
}

set rep_count [add_vdw_rep $mol_id $rep_count [$DNA_mother text] $DNA_rad 26]
set rep_count [add_vdw_rep $mol_id $rep_count [$DNA_left text] $DNA_rad 12]
set rep_count [add_vdw_rep $mol_id $rep_count [$DNA_right text] $DNA_rad 27]
set rep_count [add_vdw_rep $mol_id $rep_count [$ori text] $ori_rad 1]
set rep_count [add_vdw_rep $mol_id $rep_count [$ter text] $ter_rad 3]
set rep_count [add_vdw_rep $mol_id $rep_count [$fork text] $fork_rad 27]
set rep_count [add_vdw_rep $mol_id $rep_count [$anchor text] $anchor_rad 16]
set rep_count [add_vdw_rep $mol_id $rep_count [$hinge text] $hinge_rad 8]
set rep_count [add_vdw_rep $mol_id $rep_count [$ribo text] $ribo_rad 13]
set rep_count [add_vdw_rep $mol_id $rep_count [$bdry text] $bdry_rad 2]
mol modmaterial [expr $rep_count - 1] $mol_id Transparent

display resetview

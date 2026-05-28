#  Reproduces the look of Figure 4 of the Bentopy tutorial.
#  Usage:   vmd -e vis_t1.tcl solvated_system.gro
#  or, if VMD is already open with the file loaded:
#           source vis_t1.tcl

# Remove the default "all / Lines" representation
mol delrep 0 top

# --- Proteins: green QuickSurf surface ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 7
mol selection {not resname W and not resname NA and not resname CL}
mol material AOShiny
mol addrep top

# --- Water: translucent blue surface ---
mol representation QuickSurf 2.5 1.0 2.0 1.0
mol color ColorID 0
mol selection {resname W}
mol material Transparent
mol addrep top

# --- Periodic-boundary box ---
pbc box -color white -width 2

# --- Display settings ---
color Display Background white
display projection Orthographic
display depthcue off
display rendermode GLSL
display ambientocclusion on
display aoambient    0.80
display aodirect     0.30
display shadows      on
axes location off

# Centre the view
mol top top
display resetview

#  Reproduces Figure 5 of the Bentopy tutorial:
#     gray membrane, green lysozymes (lyz) in bulk solvent,
#     blue ubiquitins (ubq) near the membrane, translucent water.
#
#  Usage:  vmd -e vis_t2.tcl solvated_system.gro
#  or:     source vis_t2.tcl    (if VMD already has the file loaded)

mol delrep 0 top

# --- POPC membrane: gray surface ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 2
mol selection {resname POPC}
mol material AOShiny
mol addrep top

# --- Lysozyme (lyz tag): light-green surface ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 7
mol selection {resname lyz}
mol material AOShiny
mol addrep top

# --- Ubiquitin (ubq tag): blue surface ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 0
mol selection {resname ubq}
mol material AOShiny
mol addrep top

# --- Water: very transparent ---
mol representation QuickSurf 2.5 1.0 2.0 1.0
mol color ColorID 6
mol selection {resname W}
mol material Transparent
mol addrep top

pbc box -color black -width 2

color Display Background white
display projection Orthographic
display depthcue off
display rendermode GLSL
display ambientocclusion on
display aoambient    0.80
display aodirect     0.30
display shadows      on
axes location off

mol top top
display resetview
rotate x by -75   ;# look along the membrane normal so it's visible as a slab

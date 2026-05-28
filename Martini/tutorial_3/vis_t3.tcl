#  Reproduces Figure 6 of the Bentopy tutorial:
#     two POPC bilayers (gray) creating compartments A and B,
#     green lysozymes (lyz) confined to compartment A,
#     blue ubiquitins (ubq) confined to compartment B near the membrane,
#     translucent water everywhere.
#
#  Usage:  vmd -e vis_t3.tcl solvated_system.gro
#  or:     source vis_t3.tcl    (if VMD already has the file loaded)

mol delrep 0 top

# --- Double POPC membrane: gray surface ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 2
mol selection {resname POPC}
mol material AOShiny
mol addrep top

# --- Lysozyme (lyz tag, compartment A): green ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 7
mol selection {resname lyz}
mol material AOShiny
mol addrep top

# --- Ubiquitin (ubq tag, compartment B near membrane): blue ---
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol color ColorID 0
mol selection {resname ubq}
mol material AOShiny
mol addrep top

# --- Water: very transparent (skip if too slow) ---
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
rotate x by -75   ;# look along membrane normal so the two bilayers are visible

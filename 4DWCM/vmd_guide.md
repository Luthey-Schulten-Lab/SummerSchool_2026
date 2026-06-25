# Visualizing 4DWCM Trajectory

This guide explains how to visualize 4DWCM model trajectories using [VMD 2](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD).

**OOD Desktop + VirtualGL setup:** complete [vmd_guide.md](../vmd_guide.md) sections **1–2** first (login, Desktop allocation, `module load vmd`, `vglrun`).

Pre-computed spatial trajectories live in the shared bgvl data folder. Load them from the OOD Desktop terminal after VMD is running.

## 1. Launch VMD with ffmpeg (4DWCM)

If you completed [vmd_guide.md](../vmd_guide.md) §2, VirtualGL and VMD are already loaded. Add **ffmpeg** for movie export, then start VMD:

```bash
module load ffmpeg
vglrun -d egl vmd
```

## 2. Load the trajectory

Pre-computed spatial trajectories live in the shared bgvl data folder:

```
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/MinCell.lm
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/DNA/chromosome.lammpstrj
```

Load in the traj by source TCL script:

```bash
cd /projects/bgvl/SummerSchool_2026/4DWCM/render/
source load_and_sync.tcl
```

The given TCL script load in the frames at a per-minute frenquency, a downsample of the actual 4D simulation for easy visualization.

By default, all the particles in both files are rendered as Points, which are not pretty and intuitive to understand.

## 3. Set up the representations and make a movie

In 4DWCM, so many things can be visualized. Use the given TCL script

```bash
source representations.tcl
```

to set up the following representations as in the 4DWCM, Cell paper:
Left and right chromosomes as blue and red polymers of 10 bp/bead;  
Membrane as a green shell of 10 nm;
Ribosomes as spheres, yellow for actively translating and purple for free.

The membrane is drawn as a **half cut-away** (lower hemisphere) so the chromosomes and ribosomes inside stay visible.

https://github.com/user-attachments/assets/257d09f5-e3fc-4b90-82ac-8b0f173c0b21

[Download mincell.mp4](render/mincell.mp4)

Finally, let's make a cool movie to see the dynamics over the cell cycle. Personalize your `$MovieName`. 

```bash
set movie_name $MovieName
source make_movie.tcl 
```

Run this command on your laptop to download the Movie:

```bash
scp $DeltaUserName@login.delta.ncsa.illinois.edu:/projects/bgvl/SummerSchool_2026/4DWCM/render/$MovieName.mp4 ./
```

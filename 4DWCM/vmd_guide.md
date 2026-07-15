# Visualizing 4DWCM Trajectory

This guide explains how to visualize 4DWCM model trajectories using [VMD 2](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD).

VMD and the [VMD plugin](https://github.com/Luthey-Schulten-Lab/LMVMDPlugin) are pre-installed on NCSA Delta HPC.

Since VMD requires a graphical interface, we'll use Open OnDemand's Desktop interactive app to access a Linux GUI for launching VMD and viewing trajectories.

## 1. Initialize the OOD Interactive Session

1. Navigate to the [Open OnDemand dashboard](https://openondemand.delta.ncsa.illinois.edu/pun/sys/dashboard).

2. Log in through CILogon with your NCSA username, password, and Duo MFA.

3. Open the Interactive Apps menu and click **Desktop**.

4. Configure the job settings and click Launch:
   - Container image: keep default
   - Account: `bgvl-delta-gpu`
   - Partition: `GPUA100x4`
   - Duration: `00-02:00:00`
   - Reservation: leave empty if none
   - CPUs: `16`
   - RAM: `64G`
   - GPUs: `1`

5. Wait for the job status to change from "starting" to "running" in My Interactive Sessions. 

   <img src="https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-starting.png" alt="starting" width="300">

   Click "Connect to Desktop" to access the Linux graphical interface.

   <img src="https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-connect.png" alt="running" width="300">

## 2. Group leaders: stage your run (terminal, before VMD)

**Do this once in an OOD Desktop terminal before launching VMD.** Only the group leader can run this step — you must move your own `4dwcm_run`.

`Group1`, `Group2`, and `Group3` already exist under `/projects/bgvl/Groups_4DWCM/`. After your `4dwcm_7200` job finishes, **move** (do not copy) your run into your group folder. This avoids duplicating ~40 GB on `/projects/bgvl`.

```bash
GROUP=Group1   # Group1, Group2, or Group3 — use your assigned group

mv /projects/bgvl/${USER}/4dwcm_run \
   /projects/bgvl/Groups_4DWCM/${GROUP}/
```

Members then visualize from:

```
/projects/bgvl/Groups_4DWCM/Group1/4dwcm_run/Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200
```

Skip this section if you are using the shared pre-computed trajectory (`Mar31_1`) below.

## 3. Launch VMD

Open a terminal in the OOD Desktop, then launch ffmpeg (movie maker) and VMD:

```bash
source /projects/bgvl/SummerSchool_2026/DNA/files/VirtualGL/setup_env.sh
module use /projects/bgvl/alfiaparvez/modulefiles
module load ffmpeg
module load vmd/2.0.0
vglrun -d egl vmd
```

## 4. Load the trajectory

In the **VMD TkConsole**, set `indir` to your trajectory directory, then source the load script.

**Option A — shared pre-computed trajectory** (all participants, no group run needed):

```
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/MinCell.lm
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/DNA/chromosome.lammpstrj
```

```tcl
cd /projects/bgvl/SummerSchool_2026/4DWCM/render/
set indir /projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1
source load_and_sync.tcl
```

**Option B — your group's own simulation** (after the leader completed §2). Replace `Group1` with your group (`Group1`, `Group2`, or `Group3`):

```tcl
cd /projects/bgvl/SummerSchool_2026/4DWCM/render/
set indir /projects/bgvl/Groups_4DWCM/Group1/4dwcm_run/Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200
source load_and_sync.tcl
```

The script loads frames at a per-minute frequency, a downsample of the actual 4D simulation for easy visualization.

By default, all the particles in both files are rendered as Points, which are not pretty and intuitive to understand.

## 5. Set up the representations and make a movie

In 4DWCM, so many things can be visualized. In the VMD TkConsole:

```tcl
source representations.tcl
```

This sets up the following representations as in the 4DWCM, Cell paper:
Left and right chromosomes as blue and red polymers of 10 bp/bead;  
Membrane as a green shell of 10 nm;
Ribosomes as spheres, yellow for actively translating and purple for free.

The membrane is drawn as a **half cut-away** (lower hemisphere) so the chromosomes and ribosomes inside stay visible.

https://github.com/user-attachments/assets/257d09f5-e3fc-4b90-82ac-8b0f173c0b21

[Download mincell.mp4](render/mincell.mp4)

Finally, let's make a cool movie to see the dynamics over the cell cycle. Personalize your `$MovieName`, and set `outdir` to where the rendered image frames should be written (default `/projects/bgvl/SummerSchool_2026/4DWCM/render/movies/$USER`):

```tcl
set movie_name $MovieName
set outdir /u/$USER/SummerSchool_2026/4DWCM/
source make_movie.tcl
```

Run this command on your laptop to download the Movie to local. Replace `USER` and `MovieName`.

```bash
scp $DeltaUserName@login.delta.ncsa.illinois.edu:/u/$USER/SummerSchool_2026/4DWCM/$MovieName.mp4 ./
```

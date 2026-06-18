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

## 2. Launch VMD

Open a terminal, go to your cloned RDME folder, and launch ffmpeg (movie maker) and VMD:

```bash
source /projects/bgvl/SummerSchool_2026/DNA/files/VirtualGL/setup_env.sh
module load ffmpeg
module load vmd
```

## 3. Load the trajectory

Pre-computed spatial trajectories live in the shared bgvl data folder:

```
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/MinCell.lm
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/DNA/chromosome.lammpstrj
```

Load in the traj by source TCL script:

1. Open TK console in VMD by clicking Plugins-TK console.
2. Source TCL script

```bash
cd render/
source load_and_sync.tcl
```

By default, all the particles in both files are rendered as Points, which are not pretty and intuitive to understand.

## 4. Set up the representations and make a movie

In 4DWCM, so many things can be visualized. Use the given TCL script

```bash
source representations.tcl
```

to set up the following representations as in the 4DWCM, Cell paper:
Left and right chromosomes as blue and red polymers of 10 bp/bead;  
Membrane as a green shell of 10 nm;
Ribosomes as spheres, yellow for actively translating and purple for free.

The membrane is drawn as a **half cut-away** (lower hemisphere) so the chromosomes and ribosomes inside stay visible.

Finally, let's make a cool movie to see the dynamics over the cell cycle. Personalize your `$MovieName`. 

```bash
set movie_name $MovieName
source make_movie.tcl 
```

This ray-traces every frame to `render/frames/` and encodes it. For an `mp4`, run
`module load ffmpeg` **before** launching VMD; otherwise you get a `mincell.gif`.

Run this command on your laptop to download the Movie:

```bash
scp $DeltaUserName@login.delta.ncsa.illinois.edu:/projects/bgvl/SummerSchool/4DWCM/render/$MovieName ./
```
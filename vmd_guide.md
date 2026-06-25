# Visualizing trajectories with VMD on Delta

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is used across the Summer School to view **DNA/LAMMPS** and **4DWCM** trajectories. VMD needs a **graphical desktop** — it does not run in the QCB Gateway Jupyter notebooks.

Use the **Open OnDemand (OOD) Desktop** interactive app on Delta. Complete **sections 1–2 below once**, then follow the module-specific guide for your trajectory format.

| Module | Trajectory type | Continue at |
|--------|-----------------|-------------|
| **DNA** | LAMMPS `.lammpstrj` | [DNA/README.md §12](DNA/README.md#12-visualization-with-vmd) |
| **4DWCM** | Lattice Microbes `.lm` + DNA | [4DWCM/vmd_guide.md](4DWCM/vmd_guide.md) |

---

## 1. Initialize the OOD Interactive Session

1. Navigate to the [Open OnDemand dashboard](https://openondemand.delta.ncsa.illinois.edu/pun/sys/dashboard).

2. Log in through **CILogon** with your NCSA username, password, and Duo MFA.

3. Open the **Interactive Apps** menu and click **Desktop**.

4. Configure the job settings and click **Launch**:

   | Setting | Value |
   |---------|--------|
   | Container image | keep default |
   | Account | `bgvl-delta-gpu` |
   | Partition | `cpu-interactive` |
   | Duration | `00-00:30:00` |
   | Reservation | leave empty if none |
   | CPUs | `16` |
   | RAM | `64GB` |
   | GPUs | `1` |

5. Wait for the job status to change from **starting** to **running** in **My Interactive Sessions**.

   <img src="https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-starting.png" alt="starting" width="300">

   Click **Connect to Desktop** to access the Linux graphical interface.

   <img src="https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-connect.png" alt="running" width="300">

> [!NOTE]
> Request **GPUs: 1** on the Desktop form. VMD uses VirtualGL to render on the GPU; without a GPU the session falls back to software rendering (`llvmpipe`) and large trajectories are very slow to rotate.

---

## 2. Preprocess trajectory and load VMD

VMD runs on the **OOD Desktop** (not Jupyter) and is launched through **VirtualGL** so it renders on the node's GPU — the default software renderer (`llvmpipe`) makes the full ~96k-bead cell sluggish to rotate. Make sure your Desktop session requested a **GPU** (`GPUs: 1`), then set up VirtualGL (a user-space build is provided — no install or root needed) and load VMD:

```bash
source /projects/bgvl/SummerSchool_2026/DNA/files/VirtualGL/setup_env.sh
module use /projects/bgvl/alfiaparvez/modulefiles
module load vmd/2.0.0
```

Launch VMD on the GPU:

```bash
vglrun -d egl vmd
```

**Confirm the GPU is active:** VMD's startup log (or `display glinfo` in the Tk Console) should report `OpenGL renderer: NVIDIA A100...`, *not* `llvmpipe`. If `vglrun -d egl` cannot reach the GPU, list devices with `eglinfo -e` and pass one explicitly, e.g. `vglrun -d /dev/dri/card0 vmd`.

You can then visualize a **pre-run sample trajectory** (no simulation needed) or **your own run** — see the module guide:

- **DNA (LAMMPS):** [DNA/README.md §12](DNA/README.md#12-visualization-with-vmd) — sample `full_model.lammpstrj` or your `DNA_SummerSchool_2026/data/` output
- **4DWCM (RDME + DNA):** [4DWCM/vmd_guide.md](4DWCM/vmd_guide.md) — shared `MinCell.lm` trajectories and render scripts

# 4D Whole Cell (Delta Gateway) — Optimized stack

**4D Whole-Cell Model (4DWCM)** simulates the genetically minimal cell **JCVI-syn3A** in 3D over time — coupling **RDME**, **CME**, **ODE** metabolism, and **Brownian-dynamics** DNA (btree_chromo 2.0 + Kokkos LAMMPS on GPU).


---

## Which workflow should I use?

| | **JupyterHub** (Section 1) | **SSH + Slurm** (Section 2) |
|---|---------------------------|----------------------------|
| **Best for** | Short term run | Full **7200 s** cell-cycle production run |
| **Biological time** | 10 s test (users can change)  | 7200 s (~25 h wall) |
| **How you run** | Notebook cells | `sbatch` batch job |
| **Tutorial file** | [`Tutorial_4dwcm.ipynb`](Tutorial_4dwcm.ipynb) | [`4DWCM_ssh/`](4DWCM_ssh/) bash scripts |

Pick **one** path below.

---

## Shared files on Delta (read-only — copy before you run)

Simulation code and scripts live on the shared `bgvl` project at `/projects/bgvl/containers/4DWCM_Gateway/`. **Copy** into your own folder — nothing you run writes back here. Access requires the **bgvl** allocation (`delta_bgvl` group).

| Path | What it is |
|------|------------|
| `.../Optimize_4DWCM_Minimal_Cell/` | Simulation **source** — copy into your workspace or `4dwcm_run/` |
| `.../4DWCM_ssh/` | **Bash scripts** for SSH/Slurm — copy into `4dwcm_run/` |
| `.../4dcell_delta_btree2.sif` | **Container image** — used by Slurm jobs only (not copied) |

---

## 1. JupyterHub — short interactive test

### 1.1 Connect

```bash
ssh -L 8000:dt-svc-bbkw01.hsn.cm.delta.internal.ncsa.edu:8000 YOUR_NCSA_USERNAME@login.delta.ncsa.illinois.edu
```

Open https://dt-svc-bbkw01.delta.ncsa.illinois.edu:8000/hub/org/ → **QCB Gateway** → **CI Logon**.

> Gateway access: contact Alfia Parvez at alfiap@illinois.edu if not approved.

### 1.2 Allocate resources

| Setting | Value |
| --- | --- |
| **Allocation** | **A100 GPU - up to 8 (bgvl-delta-gpu)** — **Batch** |
| **CPUs** | **8** |
| **GPU Environment** | **`4DCell Optimized`** |
| **GPUs** | **2** |
| **Memory** | **64 GB**+ |
| **Time limit** | **48 hours** |
| **Kernel** | **LM 2.5 (Python 3.7)** |

### 1.3 Run the notebook

If you cloned this repository on the Gateway:

```bash
cd /home/user/workspace/SummerSchool_2026/4DWCM/4DWCM_Simulation
```

Open [`Tutorial_4dwcm.ipynb`](Tutorial_4dwcm.ipynb) and run all cells.

- **Section 1** copies source from `/projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell` into your workspace.
- **Section 3** starts the 10 s test simulation.
- **Section 4** shows how to tail the log and inspect outputs.
- **Section 5** (optional) restarts if the run crashed or checkpointed early.

### 1.4 Where JupyterHub outputs go

| What | Path |
|------|------|
| Your code copy | `/home/user/workspace/Optimize_4DWCM_Minimal_Cell/` |
| Science data | `.../Data/4dwcm_10s/` |
| Run log | `.../logs/run_4dwcm_10s.log` |

Key files in `Data/4dwcm_10s/`: `counts_and_fluxes.csv`, `DNA/`, `restart_files/`, `*.lm`.

---

## 2. SSH + Slurm — full production run

### 2.1 Log in

```bash
ssh YOUR_NCSA_USERNAME@login.delta.ncsa.illinois.edu
```

Complete DUO authentication, then:

```bash
cd /projects/bgvl/$USER
```

On Delta, `$USER` is your NCSA username (e.g. login as `jdoe` → `/projects/bgvl/jdoe`).

### 2.2 Copy source and scripts into your folder

```bash
mkdir -p 4dwcm_run/logs

cp -r /projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell 4dwcm_run/
cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_7200.sh 4dwcm_run/
cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_restart.sh 4dwcm_run/
chmod +x 4dwcm_run/launch_4dwcm_7200.sh 4dwcm_run/launch_4dwcm_restart.sh
```

Or run the helper (same result):

```bash
bash /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/copy_4dwcm_bundle.sh
```

Your folder after copying:

```
/projects/bgvl/$USER/4dwcm_run/
├── launch_4dwcm_7200.sh       ← new run
├── launch_4dwcm_restart.sh    ← resume / extend
├── logs/                      ← Slurm stdout/stderr
└── Optimize_4DWCM_Minimal_Cell/
    ├── Whole_Cell_Minimal_Cell.py
    ├── input_data/
    └── Data/                  ← created when the job runs
```

### 2.3 Submit a new run

```bash
cd /projects/bgvl/$USER/4dwcm_run
sbatch launch_4dwcm_7200.sh
```

Defaults: **7200 s** biological time, seed **13**, output folder **`4dwcm_7200`**, checkpoint at **47.5 h** wall time (`-mh`).

Custom run:

```bash
sbatch launch_4dwcm_7200.sh 7200 13 my_run_name 47.5
#                                      ↑ output folder name (Data/my_run_name/)
```

| Slurm resource | Value |
|----------------|--------|
| Account | `bgvl-delta-gpu` |
| Partition | `gpuA100x4` |
| GPUs | **2** |
| Wall time | **48:00:00** |

The job uses the shared container image at `/projects/bgvl/containers/4DWCM_Gateway/4dcell_delta_btree2.sif` (not copied).

### 2.4 Monitor

```bash
squeue -u $USER
tail -f /projects/bgvl/$USER/4dwcm_run/logs/4dwcm_7200-<JOBID>.out
```

### 2.5 Where SSH outputs go

All outputs stay under `/projects/bgvl/$USER/4dwcm_run/`:

| What | Path (default run) |
|------|---------------------|
| Slurm log | `logs/4dwcm_7200-<JOBID>.out` |
| Science data | `Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200/` |

```bash
ls /projects/bgvl/$USER/4dwcm_run/Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200/
```

| File / folder | Purpose |
|---------------|---------|
| `counts_and_fluxes.csv` | Metabolism counts and fluxes |
| `DNA/` | Chromosome / LAMMPS outputs |
| `restart_files/` | Checkpoints for resume |
| `*.lm` | RDME lattice trajectories |

### 2.6 Resume or extend a run (optional)

Use `launch_4dwcm_restart.sh` when you need **more biological time** or must **recover after a crash or checkpoint**.

`Restart_Whole_Cell_Minimal_Cell.py` takes the **same arguments** as the main script. It only works if `Data/<outputDir>/` already contains a complete prior run — use the **same** `-od` as the original (e.g. `4dwcm_7200`).

**`-t` is additional time, not the new total.** If you reached **3600 s** and want **6000 s** total, pass **2400**.

```bash
cd /projects/bgvl/$USER/4dwcm_run
sbatch launch_4dwcm_restart.sh 2400 13 4dwcm_7200 47.5
```

| Argument | Meaning |
|----------|---------|
| `2400` | Additional biological seconds |
| `13` | Same seed as the original run |
| `4dwcm_7200` | Same output folder as the original run |
| `47.5` | Wall-clock checkpoint limit (`-mh`, hours) |

Science data appends to the same `Data/4dwcm_7200/` folder. Slurm log: `logs/4dwcm_restart-<JOBID>.out`.

---

## Files in this folder

| File | Used by |
|------|---------|
| [`Tutorial_4dwcm.ipynb`](Tutorial_4dwcm.ipynb) | JupyterHub (Section 1) |
| [`4DWCM_ssh/launch_4dwcm_7200.sh`](4DWCM_ssh/launch_4dwcm_7200.sh) | SSH — new 7200 s run |
| [`4DWCM_ssh/launch_4dwcm_restart.sh`](4DWCM_ssh/launch_4dwcm_restart.sh) | SSH — resume / extend |
| [`4DWCM_ssh/copy_4dwcm_bundle.sh`](4DWCM_ssh/copy_4dwcm_bundle.sh) | SSH — optional copy helper |

---

## References

- [4D Minimal Cell (site)](https://minimalcell4d.web.illinois.edu/home/)
- [Thornburg *et al.*, *Cell* 2026](https://www.cell.com/cell/fulltext/S0092-8674(26)00174-1)
- [4DWCM_ssh README](4DWCM_ssh/README.md)
- [4DWCM analysis tutorial](../README.md)

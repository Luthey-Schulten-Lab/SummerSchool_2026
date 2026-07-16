# 4D Whole Cell (Delta) — Optimized stack

**4D Whole-Cell Model (4DWCM)** simulates the genetically minimal cell **JCVI-syn3A** in 3D over time — coupling **RDME**, **CME**, **ODE** metabolism, and **Brownian-dynamics** DNA (btree_chromo 2.0 + Kokkos LAMMPS on GPU).

Run via **SSH + Slurm** on Delta: copy simulation source and launch scripts into your project folder, then `sbatch`.

---

## 1. Log in

```bash
ssh YOUR_NCSA_USERNAME@login.delta.ncsa.illinois.edu
```

Complete DUO authentication, then:

```bash
cd /projects/bgvl/$USER
```

On Delta, `$USER` is your NCSA username (e.g. login as `jdoe` → `/projects/bgvl/jdoe`).

---

## 2. Copy source and scripts into your folder

**Simulation source** lives on the shared `bgvl` bundle. **Launch scripts** are in this repository under [`4DWCM_ssh/`](4DWCM_ssh/).

After cloning `SummerSchool_2026`:

```bash
bash SummerSchool_2026/4DWCM/4DWCM_Simulation/4DWCM_ssh/copy_4dwcm_bundle.sh
```

Or copy manually:

```bash
mkdir -p 4dwcm_run/logs

cp -r /projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell 4dwcm_run/
cp SummerSchool_2026/4DWCM/4DWCM_Simulation/4DWCM_ssh/launch_4dwcm_7200.sh 4dwcm_run/
cp SummerSchool_2026/4DWCM/4DWCM_Simulation/4DWCM_ssh/launch_4dwcm_restart.sh 4dwcm_run/
chmod +x 4dwcm_run/launch_4dwcm_7200.sh 4dwcm_run/launch_4dwcm_restart.sh
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

---

## 3. Submit a new run

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

---

## 4. Monitor

```bash
squeue -u $USER
tail -f /projects/bgvl/$USER/4dwcm_run/logs/4dwcm_7200-<JOBID>.out
```

---

## 5. Where outputs go

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

---

## 6. Resume or extend a run (optional)

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

| File | Purpose |
|------|---------|
| [`4DWCM_ssh/launch_4dwcm_7200.sh`](4DWCM_ssh/launch_4dwcm_7200.sh) | New 7200 s run |
| [`4DWCM_ssh/launch_4dwcm_restart.sh`](4DWCM_ssh/launch_4dwcm_restart.sh) | Resume / extend |
| [`4DWCM_ssh/copy_4dwcm_bundle.sh`](4DWCM_ssh/copy_4dwcm_bundle.sh) | Copy helper |
| [`4DWCM_ssh/setup_4dwcm_run.sh`](4DWCM_ssh/setup_4dwcm_run.sh) | Verify run folder |

---

## References

- [4D Minimal Cell (site)](https://minimalcell4d.web.illinois.edu/home/)
- [Thornburg *et al.*, *Cell* 2026](https://www.cell.com/cell/fulltext/S0092-8674(26)00174-1)
- [4DWCM_ssh README](4DWCM_ssh/README.md)
- [4DWCM analysis tutorial](../README.md)
- [Minecraft world generator](../minecraft/README.md)

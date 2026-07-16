# Minecraft world generator (4DWCM)

Turn a 4DWCM `MinCell.lm` trajectory into an animated Minecraft world you can explore on your laptop.

**Terminal only** — no extra setup required for participants.

## What you get

The script **always writes to your own folder** on Delta (`$USER` = whoever is logged in):

```
/projects/bgvl/$USER/minecraft/Minecraft_Cell_Animation_World/
```

Everyone reads simulation data from the **same** group folder (`Groups_4DWCM`) or demo trajectory, but each person gets their **own** copy of the Minecraft world to download.

Download your folder to your laptop and copy it into your local Minecraft `saves/` directory.

## Quick start

SSH to Delta (or open a terminal on the QCB Gateway), then:

```bash
# Where the script lives (instructor copy or your git clone)
MC=/projects/bgvl/alfiaparvez/SummerSchool_2026/4DWCM/minecraft

# Demo — shared pre-computed trajectory (no group run needed)
bash $MC/run_generate_world.sh --shared

# Group run — after your leader moved 4dwcm_run to Groups_4DWCM
bash $MC/run_generate_world.sh --group Group1
```

When finished, the script prints your output path, e.g. `/projects/bgvl/alfiaparvez/minecraft/Minecraft_Cell_Animation_World/` (your username, not anyone else's).

Full run (~30 frames) takes roughly 20–30 minutes. Quick test:

```bash
bash $MC/run_generate_world.sh --shared --frames 1
```

## Prerequisites

> **Already done:** the shared Python environment is installed at `/projects/bgvl/SummerSchool_2026/conda-envs/minecraft`. Participants can skip setup and run `run_generate_world.sh` directly.

| Step | Who | Status |
| --- | --- | --- |
| Shared Python env | Instructor | **Done** — no action needed |
| Completed 4DWCM run | Your group | Slurm job finished; `MinCell.lm` exists (or use `--shared` demo) |
| Group staging (optional) | Group leader | Move run to `Groups_4DWCM` (see below) |

Participants do **not** need to install packages or run `setup_shared_env.sh`.

## Input data — two options

### Option A — Demo (`--shared`)

Uses the shared trajectory at:

```
/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/MinCell.lm
```

Use this to try the workflow before your group run finishes.

### Option B — Group run (`--group Group1`)

After your group leader stages the simulation (see [`4DWCM/vmd_guide.md`](../vmd_guide.md) §2):

```bash
GROUP=Group1   # Group1, Group2, or Group3
mv /projects/bgvl/${USER}/4dwcm_run /projects/bgvl/Groups_4DWCM/${GROUP}/
```

The script reads:

```
/projects/bgvl/Groups_4DWCM/Group1/4dwcm_run/Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200/MinCell.lm
```

Any group member can run the generator; output still goes to **your** `/projects/bgvl/$USER/minecraft/`.

## Options

```bash
bash run_generate_world.sh --help

--shared              Use Mar31_1 demo trajectory
--group Group1        Read from Groups_4DWCM/Group1
--run-name NAME     Data subfolder (default: 4dwcm_7200)
--frames N            Animation frames (default: 30)
```

## Play in Minecraft (laptop)

1. Download `/projects/bgvl/$USER/minecraft/Minecraft_Cell_Animation_World/` to your laptop (e.g. `scp -r`).
2. Copy the folder into your Minecraft saves directory:
   - **Windows:** `%appdata%\.minecraft\saves\`
   - **macOS:** `~/Library/Application Support/minecraft/saves/`
   - **Linux:** `~/.minecraft/saves/`
3. Launch Minecraft, open the world, and use the in-game animation controls.

## Instructor setup (already complete)

The shared env was created with `setup_shared_env.sh` and lives at:

```
/projects/bgvl/SummerSchool_2026/conda-envs/minecraft   # Python 3.10
```

Re-run `setup_shared_env.sh` only if the env is missing or packages need updating.

## Troubleshooting

| Problem | Fix |
| --- | --- |
| `Shared env missing` | Env should already exist; if not, ask instructor to re-run `setup_shared_env.sh` |
| `MinCell.lm not found` | Check group staging or use `--shared` for demo |
| Wrong group path | Use `--group Group1` (or `Group2` / `Group3`) |
| Re-run after fixing input | Safe to re-run; overwrites `lm_sim*.nbt` in your output folder |

## Files in this folder

| File | Purpose |
| --- | --- |
| `run_generate_world.sh` | **Main entry point** for participants |
| `generate_world.py` | Python logic (called by the shell script) |
| `setup_shared_env.sh` | Instructor: create shared conda env |
| `Minecraft_Cell_Animation_World/` | World template (copied to your folder on first run) |

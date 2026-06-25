# 4DWCM_ssh — Delta SSH / Slurm

Copy **source** and **launch script** into your folder, then `sbatch`. Outputs stay in your folder.

## Copy

Use your **NCSA Delta username**. On Delta, `$USER` is your login name (e.g. `jdoe` → `/projects/bgvl/jdoe`).

```bash
cd /projects/bgvl/$USER
mkdir -p 4dwcm_run/logs
cp -r /projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell 4dwcm_run/
cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_7200.sh 4dwcm_run/
cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_restart.sh 4dwcm_run/
chmod +x 4dwcm_run/launch_4dwcm_7200.sh 4dwcm_run/launch_4dwcm_restart.sh
```

## Submit (new run)

```bash
cd /projects/bgvl/$USER/4dwcm_run
sbatch launch_4dwcm_7200.sh
```

## Restart or extend a run

Use `launch_4dwcm_restart.sh` when you need more biological time or to recover after a crash/checkpoint. Same arguments as the main script; `-t` is **additional** seconds (not total). Example: ran 3600 s, want 6000 s total → pass `2400`.

```bash
cd /projects/bgvl/$USER/4dwcm_run
sbatch launch_4dwcm_restart.sh 2400 13 4dwcm_7200 47.5
```

Outputs continue in the same `Data/<OUTPUT_DIR>/` folder. Slurm log: `logs/4dwcm_restart-<JOBID>.out`.

## Where outputs go

| Output | Path |
|--------|------|
| Slurm log | `/projects/bgvl/$USER/4dwcm_run/logs/4dwcm_7200-<JOBID>.out` |
| Science data | `/projects/bgvl/$USER/4dwcm_run/Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200/` |

Nothing is written to `/projects/bgvl/containers/4DWCM_Gateway/`.

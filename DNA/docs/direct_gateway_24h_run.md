# Implementation spec: run the full DNA simulation **directly inside a long (24 h) Gateway session**

> **Status:** design note / handoff spec — **not currently wired into the tutorial.**
> The shipped tutorial uses the simpler **SSH + `sbatch`** batch-job flow (see `DNA/README.md` §5).
> This document describes the *alternative* "no-SSH, run-in-the-Gateway-session" approach that was
> prototyped and validated, so it can be implemented later by handing this file to an agent.

---

## 1. Goal and when to use this

Run the **full Syn3A chromosome simulation** (loop extrusion + replication, 54,338 DNA monomers,
~14 h on one A100) **entirely inside a single QCB Delta Gateway GPU session**, with **no SSH and no
`sbatch`**. The participant requests a long (e.g. **24 h**) Gateway session, opens a notebook, and the
simulation runs in that session.

Use this approach if the priority is "never leave the browser / never touch a terminal." The tradeoff
vs. the batch-job flow is documented in §7.

---

## 2. Key facts about the environment (verified)

- The Gateway Jupyter runs **inside an Apptainer container**. That container ships a working
  `btree_chromo` at **`/Software/btree_chromo/build/apps/btree_chromo`** — the newer **4DWCM** build
  (GPU/Kokkos).
- The container has **no Slurm client** (`sbatch`/`squeue` absent) and **cannot run Apptainer** itself.
  That is why batch submission requires SSH to a login node — and why this in-session approach exists.
- The Gateway kernel runs as a **shared service account**. It **cannot write to `/projects/bgvl/<user>`**
  (confirmed `PermissionError`). All output must go to a **writable workspace in the container home**,
  e.g. `/home/user/workspace/...`.
- The shared, read-only tutorial inputs are visible at
  **`/projects/bgvl/SummerSchool_2026/DNA/files/`**.

---

## 3. The 4DWCM API differences (why we can't just run the 2025 scripts)

The Gateway's 4DWCM `btree_chromo` is **not** API-compatible with the 2025 `translocate`-based tutorial:

| 2025 command | 4DWCM status | Replacement |
| --- | --- | --- |
| `translocate:50,T` | exists but **signature changed** — now takes **0 params** | loop params loaded from a file; loop extrusion driven by `simulator_run_loops` |
| `simulator_load_loop_params:<file>` | **supported** | keep — feeds loop parameters |
| `simulator_run_loops:<#updates>,<#steps>,<thermo>,<dump>,append?,first?` | **supported** | this is how loops are driven |
| `simulator_set_DNA_model:.../LAMMPS_DNA_model/` (CPU) | **fails**: "KOKKOS package requires a Kokkos-enabled atom_style" | use `.../LAMMPS_DNA_model_kk` |
| `load_mono_quats` / `load_ribo_quats` (ellipsoids) | **fails**: "No ellipsoids allowed with this atom style" | drop them; add `switch_twisting_angles:F` and `switch_ellipsoids:F` |
| `map_replication`, `transform:m_cw…_ccw…`, `set_initial_state`/`set_final_state`, `sync_simulator_and_system` | **supported** | keep (replication) |

These were established by iteratively debugging `DNA/files/prototype_4dwcm.inp` (a 5,000-monomer toy)
until loop extrusion + replication ran cleanly in the Gateway. That prototype file is the canonical
"known-good command flow" — read it first.

Crucially, the lab already maintains a **full-scale 4DWCM directive**:
**`DNA/files/full_model.inp`** (54,338 monomers, `LAMMPS_DNA_model_kk`, `switch_ellipsoids:F`,
growing boundary via `transform:m_cw1360_ccw1360`, `map_replication`, nested `repeat` blocks with
`simulator_run_loops`). It is written with `/mnt/` paths (for the Apptainer batch path). **Reuse it as
the single source of truth** — do not rewrite the schedule.

---

## 4. Inputs (already pre-staged — no structure generation needed)

`full_model.inp` reads these from `/projects/bgvl/SummerSchool_2026/DNA/files/`:

- `x_chain_Syn3A_chromosome_init_rep00001.bin`  (initial chromosome, 54,338 monomers)
- `2500A_bdry.bin`  (boundary beads)
- `in_BD_lengths_LAMMPS_test.txt`  (BD lengths)
- `loop_params.txt`  (SMC loop parameters)

Note: `full_model.inp` uses `load_mono_coords` + `load_bdry_coords` only — **no `load_ribo_coords`** and
**no quats** (ellipsoids are off). Do not add ribosome/quaternion loads.

---

## 5. The run recipe (exact, validated)

The Apptainer batch path binds the workspace to `/mnt` and copies inputs there first. To run directly
in the Gateway, replicate that: stage inputs into a writable workspace and rewrite every `/mnt/` in
`full_model.inp` to that workspace.

```bash
# inside the Gateway GPU session (notebook %%bash cell or JupyterLab terminal)
SRC=/projects/bgvl/SummerSchool_2026/DNA/files
WS=/home/user/workspace/dna_4dwcm_full
mkdir -p "$WS"

# 1) stage the inputs full_model.inp needs
cp -n "$SRC"/in_BD_lengths_LAMMPS_test.txt          "$WS"/
cp -n "$SRC"/x_chain_Syn3A_chromosome_init_rep00001.bin "$WS"/
cp -n "$SRC"/2500A_bdry.bin                         "$WS"/
cp -n "$SRC"/loop_params.txt                        "$WS"/

# 2) build a Gateway directive: rewrite /mnt/ -> the workspace (keeps full_model.inp untouched)
sed "s|/mnt/|$WS/|g" "$SRC"/full_model.inp > "$WS"/full_model_gateway.inp

# 3) run the in-container 4DWCM btree_chromo (LD path validated on the prototype)
export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs:/Software/LAMMPS/OMP_GPU_Kokkos/lib
/Software/btree_chromo/build/apps/btree_chromo "$WS"/full_model_gateway.inp
```

Outputs land in `$WS` (trajectory `full_model.lammpstrj`, plus `full_model*` logs/data, per
`simulator_set_output_details:/mnt/,full_model` → `$WS/,full_model`).

### Recommended notebook shape (so it survives disconnects)

Launch in the **background** with `subprocess.Popen(..., start_new_session=True)` writing stdout/stderr
to `$WS/full_model_run.log`, store the PID in `$WS/run.pid`, and provide separate cells to (a) tail the
log, (b) check the PID is alive + list output, (c) `SIGTERM` the PID. (A complete version of this
notebook was prototyped and can be reconstructed from this spec; it was removed when the tutorial
switched to the batch flow — see git history for `DNA/submit_simulation.ipynb` around commit `71d92a2~1`.)

---

## 6. Gateway session allocation

When launching the Gateway server, request:

- **Allocation:** `A100 GPU - up to 8 (bgvl-delta-gpu)` (non-interactive / Batch)
- **GPU Environment:** `4DCell (LAMMPS/LM)`
- **GPUs:** 1, **CPUs:** 8, **Memory:** 64 GB
- **Time limit:** **24 hours** (the ~14 h run must fit inside the session; a non-interactive session
  keeps running even if the browser tab is closed — reopen the Gateway to check on it)

---

## 7. Caveats (must be addressed in any implementation)

1. **Session must outlive the run.** A short (default few-hour) session is killed mid-run. Hence the
   24 h request. There is **no way to detach** the run from the session without `sbatch` (which needs
   SSH) — this is the fundamental tradeoff vs. the batch flow.
2. **Queue time is not avoided.** A long GPU allocation goes through the same Slurm queue and may
   queue *longer* than a short one. For a workshop, request/launch early.
3. **Output visibility for VMD.** Output in `/home/user/workspace/...` lives **inside the Gateway
   container** and is generally **not visible from the OOD Desktop** (a separate session) used for VMD.
   The implementation must either (a) copy `full_model.lammpstrj` to a location both can read, or
   (b) document that limitation. (The batch flow avoids this because it writes to
   `/projects/bgvl/<user>/...`, which the Desktop can read.)
4. **No writing to `/projects/bgvl/<user>` from the Gateway** (shared service account). Keep all output
   under the container home workspace.
5. **Disk/quota.** The full trajectory can be multiple GB; confirm the container-home workspace has
   room.
6. **`full_model.inp` is the source of truth.** Generate `full_model_gateway.inp` at run time via `sed`;
   never hand-edit the schedule, to avoid drift.

---

## 8. Implementation checklist for the agent

- [ ] Add a notebook (e.g. `DNA/submit_simulation_gateway.ipynb`) implementing §5 (stage → `sed` →
      background launch → monitor → stop). Use the validated `LD_LIBRARY_PATH` and paths verbatim.
- [ ] In §6-style prose, instruct requesting a **24 h** GPU Gateway session.
- [ ] Handle the **VMD visibility** caveat (§7.3) explicitly — copy the trajectory somewhere the OOD
      Desktop can read, or clearly note the limitation.
- [ ] Cross-check against `DNA/files/prototype_4dwcm.inp` (known-good command flow) and
      `DNA/files/full_model.inp` (full-scale schedule). Do not reintroduce `translocate`, ellipsoids,
      quats, or the CPU `LAMMPS_DNA_model`.
- [ ] First real full-scale launch should be **validated in an actual 24 h Gateway GPU session**
      (this could not be run from the editing environment, which is not inside the Gateway container).
- [ ] Decide whether this replaces or coexists with the SSH + `sbatch` flow in `DNA/README.md`; keep the
      two consistent.

---

## 9. Reference files

- `DNA/files/full_model.inp` — full-scale 4DWCM directive (the schedule to run).
- `DNA/files/prototype_4dwcm.inp` — validated small-scale command flow + inline notes on each fix.
- `DNA/files/loop_params.txt`, `in_BD_lengths_LAMMPS_test.txt`, `2500A_bdry.bin`,
  `x_chain_Syn3A_chromosome_init_rep00001.bin` — inputs.
- `DNA/files/run_btree_chromo.sh` — the Apptainer/SIF batch script that runs `full_model.inp` via
  `/mnt` (the path this approach mirrors, minus the SIF).
- Git history: `DNA/submit_simulation.ipynb` at `71d92a2~1` contained a complete direct-run notebook.

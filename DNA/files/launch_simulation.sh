#!/bin/bash
#SBATCH --account=bgvl-delta-gpu
#SBATCH --job-name=DNA_tutorial
#SBATCH --partition=gpuA100x4
#SBATCH --time=15:00:00
#SBATCH --mem=32g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
#SBATCH --nodes=1

# Pass your workspace as the first argument. This makes the job independent of
# your username, which differs between the Gateway (e.g. alfiap-illinois) and
# your Delta SSH login. The template is copied there by the notebook / prelaunch.
# Also send the job log into that workspace so it is visible from the Gateway:
#   sbatch --output=<WORKSPACE>/DNA_tutorial.log launch_simulation.sh <WORKSPACE>
SIM_ROOT="${1:-/projects/bgvl/${USER}/DNA_SummerSchool_2026}"

apptainer run \
    --nv \
    --writable-tmpfs \
    --no-home \
    --containall \
    --bind "${SIM_ROOT}:/mnt" \
    /projects/bgvl/SummerSchool_2026/DNA/files/DNA_summer2025.sif /bin/bash -c "cd /mnt/scripts/ && bash run_sc_chain_generation.sh && python3 run_btree_chromo.py"

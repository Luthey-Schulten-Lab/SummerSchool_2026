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

# Submit from a Delta login node (where $USER is your NCSA username) with:
#   cd /projects/bgvl/$USER && sbatch launch_simulation.sh
# Output goes to /projects/bgvl/$USER/DNA_SummerSchool_2026. You may override the
# workspace by passing it as the first argument: sbatch launch_simulation.sh <DIR>
SIM_ROOT="${1:-/projects/bgvl/${USER}/DNA_SummerSchool_2026}"

# Stage the workshop template into SIM_ROOT. This runs as the submitting Delta
# user (via Slurm), so the staged files get the correct owner/group — unlike a
# copy made on the Gateway, where the kernel runs as a shared service account
# and produces files the job cannot read. Existing files are left untouched.
TEMPLATE_DIR="/projects/bgvl/SummerSchool_2026/DNA/files/DNA_SummerSchool_2026"
mkdir -p "${SIM_ROOT}"
if [[ ! -d "${SIM_ROOT}/scripts" ]]; then
    echo "Staging workshop template into ${SIM_ROOT}"
    cp -rp "${TEMPLATE_DIR}/." "${SIM_ROOT}/"
fi

apptainer run \
    --nv \
    --writable-tmpfs \
    --no-home \
    --containall \
    --bind "${SIM_ROOT}:/mnt" \
    /projects/bgvl/SummerSchool_2026/DNA/files/DNA_summer2025.sif /bin/bash -c "cd /mnt/scripts/ && bash run_sc_chain_generation.sh && python3 run_btree_chromo.py"

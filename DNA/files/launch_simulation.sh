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
#SBATCH --output=/projects/bgvl/${USER}/DNA_tutorial.log

apptainer run \
    --nv \
    --writable-tmpfs \
    --no-home \
    --containall \
    --bind /projects/bgvl/${USER}/DNA_SummerSchool_2026:/mnt \
    /projects/bgvl/SummerSchool_2026/DNA/files/DNA_summer2025.sif /bin/bash -c "cd /mnt/scripts/ && bash run_sc_chain_generation.sh && python3 run_btree_chromo.py"

#!/bin/bash
#SBATCH --account=bddt-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --time=02:00:00
#SBATCH --mem=32g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
#SBATCH --nodes=1
#SBATCH --output=/projects/bddt/${USER}/btree_chromo_workspace/full_model.log

# Set the location of the btree_chromo files
btree_chromo_files='/projects/bddt/DNA/files'

# Set the main directory for the project/workshop and the user's subdirectory
project_dir='/projects/bddt'
user_subdir=${project_dir}/${USER}

# Use the workspace
workspace_dir=${user_subdir}/btree_chromo_workspace

# Run the job inside the Apptainer container
apptainer run \
    --nv \
    --writable-tmpfs \
    --no-home \
    --containall \
    --bind ${workspace_dir}:/mnt \
    ${btree_chromo_files}/build_kokkos_image.tar_latest.sif /bin/bash -c "\
    export LD_LIBRARY_PATH=\"/usr/local/lib64:/usr/local/lib:/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs\" && \
    cd /Software/btree_chromo/build/apps && \
    ./btree_chromo /mnt/full_model.inp"


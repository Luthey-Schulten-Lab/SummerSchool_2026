#!/bin/bash

# set the location of the btree_chromo files
btree_chromo_files='/projects/bddt/DNA/files'

# set the main directory for the project/workshop and the user's subdirectory
project_dir='/projects/bddt'
user_subdir=${project_dir}/${USER}

# use the workspace
workspace_dir=${user_subdir}/btree_chromo_workspace

# use srun to launch an interactive session
# use srun to launch an interactive session
srun \
    --pty \
    --account=bddt-delta-gpu \
    --partition=gpuA100x4 \
    --time=02:00:00 \
    --mem=32g \
    --tasks-per-node=1 \
    --cpus-per-task=1 \
    --gpus-per-node=1 \
    --gpu-bind=closest \
    --nodes=1 apptainer shell \
    --nv \
    --writable-tmpfs \
    --no-home \
    --containall \
    --bind ${workspace_dir}:/mnt\
    ${btree_chromo_files}/build_kokkos_image.tar_latest.sif


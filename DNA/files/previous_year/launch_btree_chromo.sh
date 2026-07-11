#!/bin/bash

btree_chromo_files='/projects/bgvl/SummerSchool_2026/DNA/files/legacy'
project_dir='/projects/bgvl'
user_subdir=${project_dir}/${USER}
workspace_dir=${user_subdir}/btree_chromo_workspace

srun \
    --pty \
    --account=bgvl-delta-gpu \
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
    --bind ${workspace_dir}:/mnt \
    ${btree_chromo_files}/build_kokkos_image.tar_latest.sif

#!/bin/bash

btree_chromo_files='/projects/bgvl/SummerSchool_2026/DNA/files/legacy'
project_dir='/projects/bgvl'
user_subdir=${project_dir}/${USER}

workspace_dir=${user_subdir}/btree_chromo_workspace
mkdir -p "${workspace_dir}"

rsync -av --exclude='*.sif' "${btree_chromo_files}/" "${workspace_dir}/"

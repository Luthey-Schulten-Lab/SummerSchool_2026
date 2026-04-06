#! /bin/bash

# This script is used to launch VMD

cd /projects/beyi/
module load ffmpeg
module load anaconda3_cpu
source activate /projects/beyi/sw/conda/envs/vmdplugin
module load vmd/1.9.4lm
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
vmd 

#!/bin/bash
#SBATCH --job-name=WCM_4_replicates
#SBATCH --account=beyi-delta-gpu
#SBATCH --output=WCM_4replicates.out
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mem=64g
#SBATCH --partition=gpuA40x4


export HYDRA_BOOTSTRAP=fork

# Run the mpirun.sh bash file to launch WCM in parallel
apptainer exec --nv --bind /projects/beyi/$USER/ /projects/beyi/$USER/CME/summer2025.sif bash -c "source /root/miniconda3/etc/profile.d/conda.sh && conda activate lm_2.5_dev && /projects/beyi/$USER/CME/WCM/programs/mpirun.sh"
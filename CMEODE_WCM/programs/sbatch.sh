#!/bin/bash
#SBATCH --job-name=3rep_CMEODE_WCM
#SBATCH --account=bgvl-delta-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mem=32g
#SBATCH --partition=gpuA40x4
#SBATCH --mail-user=YOUR_EMAIL_ADDRESS
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=%x-%N-%j.out
#SBATCH --error=%x-%N-%j.err

export WCM=/projects/bgvl/$USER/SummerSchool_2026/CMEODE_WCM
export OUTPUT_DIR=$WCM/output_3replicates
mkdir -p "$OUTPUT_DIR/tmp"

chmod -R 777 $WCM

apptainer exec --nv \
  --bind /projects/bgvl/$USER \
  /projects/bgvl/alfiaparvez/images/4dcell_delta_btree2.sif \
  bash -c "
    source /opt/conda/etc/profile.d/conda.sh
    conda activate lm_2.5_dev
    export PYTHONPATH=/Software/Lattice_Microbes/src/pylm:\${PYTHONPATH:-}
    export HYDRA_BOOTSTRAP=fork
    # Flush each rank's stdout per line so log_<rank>.txt tracks live (the single
    # long CME run otherwise block-buffers; no effect on results).
    export PYTHONUNBUFFERED=1
    # Pin each rank to a single math thread so the 3 ranks do not oversubscribe
    # the cores via numpy/BLAS threading (each replicate is single-threaded).
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export HDF5_USE_FILE_LOCKING=FALSE
    export TMPDIR=$OUTPUT_DIR/tmp
    export PMIX_MCA_gds=hash
    export PMIX_MCA_psec=native
    cd $WCM/programs
    mpirun -np 3 python ./WCM_CMEODE_Hook.py \
      -in $WCM/input_data/ -st cme-ode -t 7200 -o 60 -hi 1 -f $OUTPUT_DIR
  "

# Make the newly generated simulation files (.lm, CSVs, logs) accessible to all.
# The chmod before the run only covers pre-existing files; this one runs AFTER the
# run so it catches everything mpirun just created. (Skipped if the job times out or
# is killed before reaching here.)
chmod -R 777 "$OUTPUT_DIR"

# ===== mpirun / python parameters =====
# for mpirun:
#   -np  number of parallel CMEODE replicates (independent cells), e.g. 3
# for python (./WCM_CMEODE_Hook.py):
#   -in  input data directory (kinetic params, SBML, genome, ...)
#   -st  simulation type, only "cme-ode" is supported
#   -t   total (biological) simulation time, integer seconds
#   -o   output (CSV flush) interval, integer seconds. The CME runs ONCE (no
#        restart); the per-hook trajectories are appended to CSV and the in-memory
#        history is trimmed every -o seconds. (-rs is a backward-compatible alias.)
#   -hi  hook interval, integer seconds (CME<->ODE communication step)
#   -f   output directory for the trajectory .csv files and log .txt files
#   Constraint: -t must be an integer multiple of -o, and -o of -hi (e.g. -t 7200 -o 60 -hi 1)

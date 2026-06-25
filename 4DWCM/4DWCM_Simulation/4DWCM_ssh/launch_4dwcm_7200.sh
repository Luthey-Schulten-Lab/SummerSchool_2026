#!/bin/bash
#SBATCH --account=bgvl-delta-gpu
#SBATCH --job-name=4dwcm_7200
#SBATCH --partition=gpuA100x4
#SBATCH --time=48:00:00
#SBATCH --mem=96g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=closest
#SBATCH --output=logs/4dwcm_7200-%j.out
#SBATCH --error=logs/4dwcm_7200-%j.err
# ---------------------------------------------------------------------------
# Optimized 4DWCM — 7200 s biological time, 2× A100, 48 h wall.
#
#   Container : /projects/bgvl/containers/4DWCM_Gateway/4dcell_delta_btree2.sif
#   Source    : $BASE/Optimize_4DWCM_Minimal_Cell (copied from containers bundle)
#
#   sbatch launch_4dwcm_7200.sh [SIM_TIME] [RNG_SEED] [OUTPUT_DIR] [MAX_HOURS]
# ---------------------------------------------------------------------------
set -euo pipefail

SIM_TIME="${1:-7200}"
RNG_SEED="${2:-13}"
OUTPUT_DIR="${3:-4dwcm_7200}"
MAX_HOURS="${4:-47.5}"

BASE="${SLURM_SUBMIT_DIR:-/projects/bgvl/${USER}/4dwcm_run}"
SIF="${SIF:-/projects/bgvl/containers/4DWCM_Gateway/4dcell_delta_btree2.sif}"
REPO="${BASE}/Optimize_4DWCM_Minimal_Cell"

mkdir -p "${BASE}/logs"

if [[ ! -f "${REPO}/Whole_Cell_Minimal_Cell.py" ]]; then
  echo "ERROR: ${REPO}/Whole_Cell_Minimal_Cell.py not found."
  echo "Copy first (from /projects/bgvl/\$USER):"
  echo "  cp -r /projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell 4dwcm_run/"
  echo "  cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_7200.sh 4dwcm_run/"
  exit 1
fi

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: container not found at ${SIF}"
  exit 1
fi

echo "############################################################"
echo "[4dwcm_7200] start    : $(date)"
echo "[4dwcm_7200] node     : $(hostname)"
echo "[4dwcm_7200] image    : ${SIF}"
echo "[4dwcm_7200] repo     : ${REPO}"
echo "[4dwcm_7200] sim/seed : ${SIM_TIME}s / ${RNG_SEED}"
echo "[4dwcm_7200] maxHours : ${MAX_HOURS}"
echo "[4dwcm_7200] GPUs     : 2 (RDME/CME GPU 0, DNA GPU 1)"
echo "[4dwcm_7200] outdir   : ${REPO}/Data/${OUTPUT_DIR}"
echo "[4dwcm_7200] slurm log: ${BASE}/logs/4dwcm_7200-${SLURM_JOB_ID:-local}.out"
echo "############################################################"

rm -rf "${REPO}/pyxbld" "${REPO}/cythonCompiledFunctions.pyx" \
       "${REPO}/cythonCompiledFunctions.c" "${REPO}/setup_tmp.py" \
       "${REPO}"/cythonCompiledFunctions*.so 2>/dev/null || true

apptainer exec \
    --nv \
    --writable-tmpfs \
    --no-home \
    --containall \
    --pwd /mnt \
    --bind "${REPO}:/mnt" \
    "${SIF}" /bin/bash -c "\
        export HOME=/tmp && \
        export TMPDIR=/tmp && \
        export XDG_CACHE_HOME=/tmp/.cache && \
        export PYTHONPYCACHEPREFIX=/tmp/.pycache && \
        export HDF5_USE_FILE_LOCKING=FALSE && \
        export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:\${LD_LIBRARY_PATH:-} && \
        export PYTHONPATH=/Software/Lattice_Microbes_2.6/src/pylm && \
        export PATH=/opt/conda/envs/lm_2.5_dev/bin:\${PATH} && \
        export DNA_GPU_ID=1 && \
        cd /mnt && \
        /opt/conda/envs/lm_2.5_dev/bin/python -u Whole_Cell_Minimal_Cell.py \
            -od ${OUTPUT_DIR} \
            -t ${SIM_TIME} \
            -cd 0 \
            -drs ${RNG_SEED} \
            -mh ${MAX_HOURS} \
            -dsd /Software/opt/"

echo "############################################################"
echo "[4dwcm_7200] DONE: $(date)"
echo "############################################################"

#!/bin/bash
#SBATCH --account=bgvl-delta-gpu
#SBATCH --job-name=4dwcm_restart
#SBATCH --partition=gpuA100x4
#SBATCH --time=48:00:00
#SBATCH --mem=96g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=closest
#SBATCH --output=logs/4dwcm_restart-%j.out
#SBATCH --error=logs/4dwcm_restart-%j.err
# ---------------------------------------------------------------------------
# Optimized 4DWCM — restart / extend an existing run (Restart_Whole_Cell_Minimal_Cell.py)
#
#   Container : /projects/bgvl/containers/4DWCM_Gateway/4dcell_delta_btree2.sif
#   Source    : $BASE/Optimize_4DWCM_Minimal_Cell (copied from containers bundle)
#
#   sbatch launch_4dwcm_restart.sh ADDITIONAL_SIM_TIME [RNG_SEED] [OUTPUT_DIR] [MAX_HOURS]
#
#   ADDITIONAL_SIM_TIME = more biological seconds to run (not total simulated time).
#   OUTPUT_DIR must match the -od used in the original run (e.g. 4dwcm_7200).
# ---------------------------------------------------------------------------
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch launch_4dwcm_restart.sh ADDITIONAL_SIM_TIME [RNG_SEED] [OUTPUT_DIR] [MAX_HOURS]"
  echo ""
  echo "  ADDITIONAL_SIM_TIME  More biological seconds to run (not total)."
  echo "                       Example: ran 3600 s, want 6000 s total → pass 2400."
  echo "  RNG_SEED             Same seed as the original run (default: 13)."
  echo "  OUTPUT_DIR           Same -od as the original run (default: 4dwcm_7200)."
  echo "  MAX_HOURS            Wall-clock checkpoint limit in hours (default: 47.5)."
  exit 1
fi

ADDITIONAL_SIM_TIME="${1}"
RNG_SEED="${2:-13}"
OUTPUT_DIR="${3:-4dwcm_7200}"
MAX_HOURS="${4:-47.5}"

BASE="${SLURM_SUBMIT_DIR:-/projects/bgvl/${USER}/4dwcm_run}"
SIF="${SIF:-/projects/bgvl/containers/4DWCM_Gateway/4dcell_delta_btree2.sif}"
REPO="${BASE}/Optimize_4DWCM_Minimal_Cell"
DATA_DIR="${REPO}/Data/${OUTPUT_DIR}"

mkdir -p "${BASE}/logs"

if [[ ! -f "${REPO}/Restart_Whole_Cell_Minimal_Cell.py" ]]; then
  echo "ERROR: ${REPO}/Restart_Whole_Cell_Minimal_Cell.py not found."
  echo "Copy first (from /projects/bgvl/\$USER):"
  echo "  cp -r /projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell 4dwcm_run/"
  exit 1
fi

if [[ ! -f "${DATA_DIR}/sim_properties.pkl" ]]; then
  echo "ERROR: no prior run found at ${DATA_DIR}/"
  echo "Restart requires Data/${OUTPUT_DIR}/ from a previous Whole_Cell_Minimal_Cell.py run"
  echo "(sim_properties.pkl, restart_files/, MinCell_restart_*.lm, etc.)."
  exit 1
fi

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: container not found at ${SIF}"
  exit 1
fi

echo "############################################################"
echo "[4dwcm_restart] start    : $(date)"
echo "[4dwcm_restart] node     : $(hostname)"
echo "[4dwcm_restart] image    : ${SIF}"
echo "[4dwcm_restart] repo     : ${REPO}"
echo "[4dwcm_restart] add time : +${ADDITIONAL_SIM_TIME}s biological"
echo "[4dwcm_restart] seed     : ${RNG_SEED}"
echo "[4dwcm_restart] maxHours : ${MAX_HOURS}"
echo "[4dwcm_restart] GPUs     : 2 (RDME/CME GPU 0, DNA GPU 1)"
echo "[4dwcm_restart] outdir   : ${DATA_DIR}"
echo "[4dwcm_restart] slurm log: ${BASE}/logs/4dwcm_restart-${SLURM_JOB_ID:-local}.out"
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
        /opt/conda/envs/lm_2.5_dev/bin/python -u Restart_Whole_Cell_Minimal_Cell.py \
            -od ${OUTPUT_DIR} \
            -t ${ADDITIONAL_SIM_TIME} \
            -cd 0 \
            -drs ${RNG_SEED} \
            -mh ${MAX_HOURS} \
            -dsd /Software/opt/"

echo "############################################################"
echo "[4dwcm_restart] DONE: $(date)"
echo "############################################################"

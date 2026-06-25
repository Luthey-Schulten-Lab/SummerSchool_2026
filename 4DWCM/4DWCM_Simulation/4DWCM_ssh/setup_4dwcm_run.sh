#!/usr/bin/env bash
# Verify /projects/bgvl/$USER/4dwcm_run is ready (run after copy).
set -euo pipefail

RUN_DIR="/projects/bgvl/${USER}/4dwcm_run"
REPO="${RUN_DIR}/Optimize_4DWCM_Minimal_Cell"
LAUNCH="${RUN_DIR}/launch_4dwcm_7200.sh"
RESTART="${RUN_DIR}/launch_4dwcm_restart.sh"

if [[ ! -f "${REPO}/Whole_Cell_Minimal_Cell.py" ]]; then
  echo "ERROR: source missing at ${REPO}"
  echo "Run copy first (from /projects/bgvl/\$USER):"
  echo "  cp -r /projects/bgvl/containers/4DWCM_Gateway/Optimize_4DWCM_Minimal_Cell 4dwcm_run/"
  echo "  cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_7200.sh 4dwcm_run/"
  echo "  cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_restart.sh 4dwcm_run/"
  exit 1
fi
if [[ ! -x "${LAUNCH}" ]]; then
  echo "ERROR: launch script missing at ${LAUNCH}"
  echo "  cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_7200.sh 4dwcm_run/"
  exit 1
fi
if [[ ! -x "${RESTART}" ]]; then
  echo "ERROR: restart script missing at ${RESTART}"
  echo "  cp /projects/bgvl/containers/4DWCM_Gateway/4DWCM_ssh/launch_4dwcm_restart.sh 4dwcm_run/"
  exit 1
fi

mkdir -p "${RUN_DIR}/logs"
echo "OK: ${RUN_DIR} is ready."
echo "New run:    cd ${RUN_DIR} && sbatch launch_4dwcm_7200.sh"
echo "Restart:    cd ${RUN_DIR} && sbatch launch_4dwcm_restart.sh <ADDITIONAL_SECS> [SEED] [OUTPUT_DIR] [MAX_HOURS]"

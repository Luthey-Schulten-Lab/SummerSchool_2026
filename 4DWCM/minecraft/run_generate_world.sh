#!/bin/bash
# Generate a Minecraft world in /projects/bgvl/$USER/minecraft/ (terminal only, no Jupyter).
#
# Demo (shared trajectory):
#   bash run_generate_world.sh --shared
#
# Group simulation (after leader moved 4dwcm_run to Groups_4DWCM):
#   bash run_generate_world.sh --group Group1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="/projects/bgvl/SummerSchool_2026/conda-envs/minecraft"
USER_WORK="/projects/bgvl/${USER}/minecraft"
WORLD_DIR="${USER_WORK}/Minecraft_Cell_Animation_World"
TEMPLATE="${SCRIPT_DIR}/Minecraft_Cell_Animation_World"
SCHEM_DIR="${WORLD_DIR}/datapacks/cellgen/data/cellgen/structure/animate"

GROUP=""
RUN_NAME="4dwcm_7200"
FRAME_NUM=30

usage() {
  sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --shared) GROUP="shared"; shift ;;
    --group) GROUP="$2"; shift 2 ;;
    --run-name) RUN_NAME="$2"; shift 2 ;;
    --frames) FRAME_NUM="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    *) echo "Unknown option: $1" >&2; usage 1 ;;
  esac
done

if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  echo "ERROR: Shared env missing. Ask instructor to run setup_shared_env.sh" >&2
  exit 1
fi

mkdir -p "${USER_WORK}"

if [[ ! -d "${WORLD_DIR}/datapacks" ]]; then
  echo "Copying world template to ${WORLD_DIR} (once)..."
  rsync -a --exclude='.git' "${TEMPLATE}/" "${WORLD_DIR}/"
fi

mkdir -p "${SCHEM_DIR}"

PY_ARGS=(--output-path "${SCHEM_DIR}" --frames "${FRAME_NUM}" --run-name "${RUN_NAME}")
if [[ -n "${GROUP}" ]]; then
  PY_ARGS+=(--group "${GROUP}")
else
  PY_ARGS+=(--user "${USER}")
fi

echo "Output directory: ${WORLD_DIR}"
"${ENV_PREFIX}/bin/python" "${SCRIPT_DIR}/generate_world.py" "${PY_ARGS[@]}"

echo ""
echo "World ready: ${WORLD_DIR}"
echo "Download that folder to your laptop and copy into Minecraft saves/."

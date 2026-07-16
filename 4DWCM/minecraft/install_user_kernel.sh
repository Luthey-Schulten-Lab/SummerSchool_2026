#!/bin/bash
# Register the shared Minecraft kernel in your OOD Jupyter (run once per user).
#
#   bash install_user_kernel.sh

set -euo pipefail

ENV_PREFIX="/projects/bgvl/SummerSchool_2026/conda-envs/minecraft"
KERNEL_NAME="bgvl-minecraft"
DISPLAY_NAME="Python (bgvl Minecraft)"

if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  echo "ERROR: Shared env not found. Ask instructor to run setup_shared_env.sh first." >&2
  exit 1
fi

"${ENV_PREFIX}/bin/python" -m ipykernel install \
  --user \
  --name="${KERNEL_NAME}" \
  --display-name="${DISPLAY_NAME}"

echo "Kernel installed: ${DISPLAY_NAME}"
echo "In Jupyter: Kernel → Change Kernel → ${DISPLAY_NAME}"

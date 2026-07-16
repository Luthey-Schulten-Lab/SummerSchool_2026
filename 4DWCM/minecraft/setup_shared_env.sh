#!/bin/bash
# Create shared conda env for generate_world.ipynb (run once from a login node).
#
#   bash setup_shared_env.sh
#
# Installs to: /projects/bgvl/SummerSchool_2026/conda-envs/minecraft
# Group: delta_bgvl (readable/executable by all bgvl users)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="/projects/bgvl/SummerSchool_2026/conda-envs/minecraft"

module load anaconda3 2>/dev/null || true
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Run: module load anaconda3" >&2
  exit 1
fi

echo "Creating shared env at: ${ENV_PREFIX}"
mkdir -p "$(dirname "${ENV_PREFIX}")"

if [[ -x "${ENV_PREFIX}/bin/python" ]]; then
  echo "Env already exists — updating packages..."
else
  conda env create -p "${ENV_PREFIX}" -f "${SCRIPT_DIR}/environment.yml" -y
fi

# Install mcschematic-plus from this repo (editable not required for tutorial)
"${ENV_PREFIX}/bin/pip" install --upgrade pip
"${ENV_PREFIX}/bin/pip" install "${SCRIPT_DIR}"

# Register Jupyter kernel inside the shared env (visible when JUPYTER_PATH is set)
"${ENV_PREFIX}/bin/python" -m ipykernel install \
  --prefix="${ENV_PREFIX}" \
  --name=bgvl-minecraft \
  --display-name="Python (bgvl Minecraft)"

# delta_bgvl group access
chgrp -R delta_bgvl "${ENV_PREFIX}" 2>/dev/null || true
chmod -R g+rX "${ENV_PREFIX}" 2>/dev/null || true
find "${ENV_PREFIX}" -type d -exec chmod g+rX {} + 2>/dev/null || true

echo ""
echo "Done. Shared env: ${ENV_PREFIX}"
echo "Python: $(${ENV_PREFIX}/bin/python --version)"
echo ""
echo "Participants — register kernel for OOD Jupyter (once per user):"
echo "  bash ${SCRIPT_DIR}/install_user_kernel.sh"
echo ""
echo "Or run the notebook without Jupyter:"
echo "  ${ENV_PREFIX}/bin/python your_script.py"

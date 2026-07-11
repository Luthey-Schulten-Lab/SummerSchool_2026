#!/bin/bash
# Build btree_chromo inside DNA_summer2025.sif.
#
# Prerequisites (one-time on bgvl): files/btree_chromo_gpu/ and files/DNA_summer2025.sif
#
# Run on a Delta login node:
#   bash /projects/bgvl/SummerSchool_2026/DNA/files/build_btree_chromo.sh
set -euo pipefail

FILES_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF="${FILES_ROOT}/DNA_summer2025.sif"
SRC="${FILES_ROOT}/btree_chromo_gpu"
BTREE_BIN="${SRC}/build/apps/btree_chromo"

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: Apptainer image not found: ${SIF}" >&2
  exit 1
fi

if [[ ! -d "${SRC}/src" ]]; then
  echo "ERROR: btree_chromo_gpu source not found: ${SRC}" >&2
  exit 1
fi

echo "FILES_ROOT=${FILES_ROOT}"
echo "Building btree_chromo from ${SRC} ..."
apptainer exec --nv --bind "${SRC}:/src" "${SIF}" /bin/bash -c '
  set -euo pipefail
  cd /src
  make clean 2>/dev/null || true
  make -j"$(nproc)"
  ls -lh build/apps/btree_chromo
  file build/apps/btree_chromo
'

echo "Build complete: ${BTREE_BIN}"

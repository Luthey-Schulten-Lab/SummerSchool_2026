#!/bin/bash
# Copy the DNA workshop launch script and template into your personal bgvl folder.
#
# Run this on a Delta LOGIN NODE (after `ssh USERNAME@login.delta.ncsa.illinois.edu`),
# where $USER is your NCSA username and /projects/bgvl/$USER is your personal folder.
#
#   bash /projects/bgvl/SummerSchool_2026/DNA/files/prelaunch_btree_chromo.sh
#
# Then:  cd /projects/bgvl/$USER && sbatch launch_simulation.sh

set -euo pipefail

SRC="/projects/bgvl/SummerSchool_2026/DNA/files"
DEST="/projects/bgvl/${USER}"

mkdir -p "${DEST}"
cp -p  "${SRC}/launch_simulation.sh"   "${DEST}/"
cp -rp "${SRC}/DNA_SummerSchool_2026"  "${DEST}/"

echo "Copied launch_simulation.sh and DNA_SummerSchool_2026/ into ${DEST}"
echo "Next:  cd ${DEST} && sbatch launch_simulation.sh"

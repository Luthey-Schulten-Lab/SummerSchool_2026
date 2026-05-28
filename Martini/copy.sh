#!/bin/bash
#
# Copy the Martini / Bentopy tutorial materials from the workshop master
# directory into your own user space on Delta:
#
#     /projects/bgvl/$USER/SummerSchool_2026/Martini/
#
# After running, open Jupyter Notebook (see ../README.md and the section
# "1. Set up the tutorial on Delta" of this README), navigate to one of:
#
#     Martini/tutorial_1/   - basic packing       (lysozyme in water)
#     Martini/tutorial_2/   - membrane packing    (Aquaporin Z + lipid bilayer)
#     Martini/tutorial_3/   - multi-compartment   (POPC vesicle)
#
# and open the .ipynb that lives in that folder.

set -euo pipefail

MASTER="/projects/bgvl/alfiaparvez/SummerSchool_2026/Martini"
DEST_PARENT="/projects/bgvl/${USER}/SummerSchool_2026"
DEST="${DEST_PARENT}/Martini"

mkdir -p "${DEST_PARENT}"

if [[ -e "${DEST}" ]]; then
    echo "[copy.sh] ${DEST} already exists - filling in any missing files only."
    cp -rn "${MASTER}/." "${DEST}/"
else
    echo "[copy.sh] copying ${MASTER}  ->  ${DEST}"
    cp -r "${MASTER}" "${DEST_PARENT}/"
fi

echo
echo "[copy.sh] done."
echo "    Open Jupyter and run:"
echo "        ${DEST}/tutorial_1/tutorial_1_basic_packing.ipynb"
echo "        ${DEST}/tutorial_2/tutorial_2_membrane_packing.ipynb"
echo "        ${DEST}/tutorial_3/tutorial_3_multi_compartment.ipynb"

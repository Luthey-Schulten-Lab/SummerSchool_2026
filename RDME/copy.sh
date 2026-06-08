#!/bin/bash
#
# RDME pre-computed data lives on bgvl (not in the git repository):
#
#   /projects/bgvl/SummerSchool_2026/RDME/data/        # counts_and_fluxes.*.csv
#   /projects/bgvl/SummerSchool_2026/RDME/trajectory/  # MinCell_*.lm for VMD
#
# Clone the repository for notebooks, then read data from those paths.
# In Jupyter:
#   !git clone https://github.com/Luthey-Schulten-Lab/SummerSchool_2026.git
#   %cd SummerSchool_2026/RDME

echo "[copy.sh] No local copy needed." >&2
echo "Clone the repo for notebooks; read data from /projects/bgvl/SummerSchool_2026/RDME/data/" >&2
exit 0

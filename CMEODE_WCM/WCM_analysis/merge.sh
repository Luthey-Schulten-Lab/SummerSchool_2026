source ~/anaconda3/etc/profile.d/conda.sh
conda activate rapids-24.04   # any env with pandas + numpy works (no cudf / GPU required)

# base directory that contains the per-launch output folders
MAIN_FOLDER="/home/enguang/Documents/Workshops/2026-STC_QCB-SummerSchool/accelerate_CMEODE/WCM_w2DCplx"

# one folder per mpirun launch, each holding counts_<n>/SA_<n>/Flux_<n>/log_<n>
OUTPUT_FOLDERS=(output_consistency output_consistency2)

python -m cudf.pandas merge_categorize.py \
    --main_folder "$MAIN_FOLDER" \
    --output_folders "${OUTPUT_FOLDERS[@]}"

# merge_categorize.py uses plain `import pandas as pd` (no cudf in the code).
# On a RAPIDS/GPU machine you can GPU-accelerate it transparently by prepending
# the accelerator -- never required, gives identical results:
#   python -m cudf.pandas merge_categorize.py --main_folder ... --output_folders ...

# merge_categorize.py
#   --main_folder     base directory containing the output folders
#   --output_folders  one or more result folders (one per mpirun launch), relative
#                     to main_folder; reps/timepoints are auto-detected
#
# Combined trajectories are written, renumbered per category, to:
#   $MAIN_FOLDER/trajs_healthy/  $MAIN_FOLDER/trajs_other/  $MAIN_FOLDER/trajs_short/
# as counts_SA_fluxes_<n>.csv (+ a copy of each log_<n>.txt).
#
# Classification metric is set by METRIC in merge_categorize.py:
#   0: full length, ATP & GTP accumulation == 0, pep > 1
#   1: full length, ATP & GTP shortage    == 0, pep > 1   (default)
#   2: ATP accumulation >= 0 up to the surface-area doubling time
#   3: pep > 0 up to the surface-area doubling time

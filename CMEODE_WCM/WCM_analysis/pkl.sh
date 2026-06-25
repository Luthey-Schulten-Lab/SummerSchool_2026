source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

# Pin math libraries to one thread so the merge_trajs worker PROCESSES do not
# oversubscribe the cores via BLAS threading.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# base directory (same one passed to merge.sh)
MAIN_FOLDER="/home/enguang/Documents/Workshops/2026-STC_QCB-SummerSchool/accelerate_CMEODE/WCM_w2DCplx"

# merged category folder(s) under MAIN_FOLDER to serialize
OUTPUT_FOLDERS=(trajs_healthy)

python pkl.py \
    --main_folder "$MAIN_FOLDER" \
    --output_folders "${OUTPUT_FOLDERS[@]}"

# pkl.py
#   --main_folder     base directory (same one passed to merge.sh)
#   --output_folders  merged category folder(s) under main_folder to serialize
#                     (default: healthy); replicate count is auto-detected
#   --n_cpus          worker processes for the parallel CSV read (default: 8)
#
# Writes WCMensemble_<main>.pkl + _x.npy/_fx.npy into $MAIN_FOLDER/<first folder>/.

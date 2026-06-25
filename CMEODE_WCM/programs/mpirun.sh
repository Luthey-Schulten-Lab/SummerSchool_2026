# The bash file to launch parallel CMEODE simulations
# Each CMEODE simulation is independent with each other, i.e. do not communicate with each other

# Activate conda enviroment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate LM_Cell

# Base directory without slash between SIM_NAME and SIM_YEAR
BASE_DIR=".."
OUTPUT_DIR="$BASE_DIR/output_consistency2/"
INPUT_DIR="$BASE_DIR/input_data/"

# Create Output Folder
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Run Simulation: 25 independent replicates (single CME run each, no restart)
mpirun -np 25 python ./WCM_CMEODE_Hook.py \
    -in "$INPUT_DIR" \
    -st cme-ode \
    -t 7200 \
    -o 60 \
    -hi 1 \
    -f "$OUTPUT_DIR"

# Input Arguments
# for mpirun:

    # -np number of parallel CMEODE replicates (independent cells), integer from 1 to nmax (here 25)

# for python:
    # -in input directory

    # -st simulation type, only support "cme-ode"

    # -t simulation time, integer numbers, in seconds

    # -o output (CSV flush) interval, integer numbers, in seconds.
    #    The CME no longer restarts: it runs once for the whole -t, and the per-hook
    #    trajectories are appended to CSV (and the in-memory history trimmed) every -o seconds.
    #    (-rs is still accepted as a backward-compatible alias of -o.)

    # -hi hook interval, integer numbers, in seconds

    # -f directory to store output trajectory .csv files and log .txt files, strings, created automatically

    # For the times, the former should be the integer multiple of the latter, e.g. -t 7200 -o 60 -hi 1

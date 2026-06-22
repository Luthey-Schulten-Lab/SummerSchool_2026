# The bash file to launch parallel CMEODE simulations
# Each CMEODE simulation is independent with each other, i.e. do not communicate with each other

# Create Output Folder
OUTPUT_DIR='../output_4replicates'
mkdir -p "$OUTPUT_DIR"

# Run Simulation
mpirun -np 4 python ./WCM_CMEODE_Hook.py -st cme-ode -t 60 -rs 60 -hi 1 -f "$OUTPUT_DIR"

# Input Arguments
# for mpirun:

    # -np numbers of parallel CMEODE simulations, integer number from 1 to nmax

# for python:

    # -st simulation type, only support "cme-ode"

    # -t simulation time, integer numbers, in seconds

    # -rs restart interval, integer numbers, in seconds

    # -hi hook interval, integer numbers, in seconds
    
    # -f directory to store output trajectory .csv files and log .txt files, strings, created automatically

    # For the times, the former should be the integer multiples of the latter e.g. -t 120 -rs 60 -wi 2 -hi 1
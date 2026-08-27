#!/bin/bash
#SBATCH --job-name=HERT_sim     # job name
#SBATCH --nodes=1               # node(s) required for the job
#SBATCH --ntasks=45             # number of tasks (cores) across all nodes
#SBATCH --partition=mak0037_std # name of partition to submit job
#SBATCH --output=test-%j.out    # output file where %j is the job ID
#SBATCH --error=test-%j.err     # error file where %j is the job ID
#SBATCH --mail-type=ALL         # will send email for begin, end, and failure
#SBATCH --mail-user=wzt0020@auburn.edu

module load geant4/11.1.2
module load cmake

# ==============================================================================
# Simulation Scaling & Energy Bins
# ==============================================================================
NUM_CORES=45
TOTAL_RUNS=9000
CURRENT_RUN=1

# Format: "START_RUN END_RUN ENERGY_MIN ENERGY_MAX"
ENERGY_BINS=(
    "1 10000 10 100"
    "10001 11000 100 200"
    "11001 12000 200 300"
    "12001 13000 300 400"
    "13001 15000 400 1000"
    "15001 17000 1000 2000"
)

# ==============================================================================
# Path Initialization 
# ==============================================================================
BASE_DIR=$(pwd)
MACRO_DIR=$(realpath ../)
RESULTS_DIR=$(realpath ../../FS_Results)
SEEDS_DIR=$(realpath ../../proton_seeds)

mkdir -p "${RESULTS_DIR}"
mkdir -p "${SEEDS_DIR}"

echo "Launching ${NUM_CORES} parallel Geant4 instances with dynamic energy bins..."

# Calculate baseline runs per core (Integer division + remainder distribution)
BASE_RPC=$(( TOTAL_RUNS / NUM_CORES ))
REMAINDER=$(( TOTAL_RUNS % NUM_CORES ))

# ==============================================================================
# Parallel Execution Loop
# ==============================================================================
for (( i=1; i<=NUM_CORES; i++ ))
do
    # 1. Determine the exact start and end run for this specific core
    CORE_START=$CURRENT_RUN
    RPC=$BASE_RPC
    if [ $i -le $REMAINDER ]; then
        RPC=$(( RPC + 1 )) # Distribute remainder evenly among early cores
    fi
    CORE_END=$(( CURRENT_RUN + RPC - 1 ))
    CURRENT_RUN=$(( CORE_END + 1 ))

    # 2. Create and enter an isolated worker directory
    WORK_DIR="${BASE_DIR}/worker_${i}"
    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}" || exit

    # 3. Construct the MultipleRun_worker.mac file dynamically
    MACRO_CONTENT="/control/shell rm -f *.txt\n"

    for BIN in "${ENERGY_BINS[@]}"; do
        # Extract values from the array string
        read BIN_START BIN_END E_MIN E_MAX <<< "$BIN"

        # Check for intersection between the Core's assignment and the Energy Bin
        OVERLAP_START=$(( CORE_START > BIN_START ? CORE_START : BIN_START ))
        OVERLAP_END=$(( CORE_END < BIN_END ? CORE_END : BIN_END ))

        if [ $OVERLAP_START -le $OVERLAP_END ]; then
            # An overlap exists! Write a loop block for this energy segment
            MACRO_CONTENT+="\n# Assigned Energy Range: ${E_MIN} to ${E_MAX} MeV\n"
            MACRO_CONTENT+="/control/alias energyMin ${E_MIN}\n"
            MACRO_CONTENT+="/control/alias energyMax ${E_MAX}\n"
            MACRO_CONTENT+="/control/alias startNumber ${OVERLAP_START}\n"
            MACRO_CONTENT+="/control/alias endNumber ${OVERLAP_END}\n"
            MACRO_CONTENT+="/control/alias stepNumber 1\n"
            MACRO_CONTENT+="/control/loop ${BASE_DIR}/SingleRunFS_worker_${i}.mac runNumber {startNumber} {endNumber} {stepNumber}\n"
        fi
    done
    MACRO_CONTENT+="exit\n"

    # Write the constructed string into the file
    echo -e "$MACRO_CONTENT" > MultipleRun_worker.mac

    # 4. Duplicate SingleRunFS.mac and strip out the hardcoded energies
    # so it inherits the energyMin/energyMax set in the MultipleRun file above.
    sed -e "s|../../proton_seeds|${SEEDS_DIR}|g" \
        -e "s|../../FS_Results|${RESULTS_DIR}|g" \
        -e "s|/control/alias energyMin|#/control/alias energyMin|g" \
        -e "s|/control/alias energyMax|#/control/alias energyMax|g" \
        "${MACRO_DIR}/SingleRunFS.mac" > "${BASE_DIR}/SingleRunFS_worker_${i}.mac"

    # 5. Execute HERT in the background
    ${BASE_DIR}/HERT MultipleRun_worker.mac > "sim_core_${i}.log" 2>&1 &
    
    echo "  -> Core ${i} launched (Runs ${CORE_START} to ${CORE_END})"
    
    # Return to base directory
    cd "${BASE_DIR}" || exit
done

# ==============================================================================
# Synchronization and Cleanup
# ==============================================================================
echo "All ${NUM_CORES} instances submitted. Waiting for completion..."
wait 

echo "Simulations complete! Cleaning up worker directories..."

for (( i=1; i<=NUM_CORES; i++ ))
do
    rm -rf "${BASE_DIR}/worker_${i}"
    rm -f "${BASE_DIR}/SingleRunFS_worker_${i}.mac"
done

echo "Done! Data safely moved to ${RESULTS_DIR}"
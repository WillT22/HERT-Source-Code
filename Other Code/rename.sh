#!/bin/bash

# Define the boundaries of the INCORRECTLY numbered files
START_RUN=12001
END_RUN=12200
OFFSET=1000
BASE_PATTERN="HERT_CADoutput_proton_1000000_Run" # Base part of the filename

# Loop through all the misnamed run numbers from 12001 to 12200
for old_run_num in $(seq $START_RUN $END_RUN); do

    # Construct the full OLD filename with the .txt extension
    old_file="${BASE_PATTERN}${old_run_num}.txt"

    # 1. Check if the file actually exists
    if [[ -f "$old_file" ]]; then

        # 2. Calculate the new, correct run number (e.g., 12001 - 1000 = 11001)
        new_run_num=$((old_run_num - OFFSET))

        # 3. Construct the full NEW filename
        new_file="${BASE_PATTERN}${new_run_num}.txt"

        # 4. Perform the rename operation
        echo "Renaming: $old_file  -->  $new_file"
        mv "$old_file" "$new_file"
    fi
done

echo "Renaming batch complete."
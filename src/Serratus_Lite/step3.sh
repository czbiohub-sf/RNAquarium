#!/bin/bash

# Step 3
# WAIT FOR HPC TO FINISH STEP 2
WORKING_DIR=$(realpath .)
for VIRUS_DIR in results/* ; do
    cd "$(realpath "${VIRUS_DIR}")"
    python3 "${WORKING_DIR}/scripts/aggregate_summarized_files.py" --dir .
done

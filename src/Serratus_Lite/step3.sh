#!/bin/bash

# Step 3
# WAIT FOR HPC TO FINISH STEP 2
for VIRUS_DIR in results/* ; do
    cd "$(realpath "${VIRUS_DIR}")"
    python3 aggregate_summarized_files.py --dir .
done

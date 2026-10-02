#!/bin/bash
REPO=$(realpath "$(dirname "$0")/../..")
export PYTHONPATH=$REPO

python -m running_from_file.run_on_directory \
    --input-dir  "$REPO/data/gib_files_filtered_oct1_2026_new_format" \
    --output-dir "$REPO/Negatives/count_GIAB_mutation_rates/results/with_fisher_fixed" \
    --workers 6
#!/bin/bash
REPO=$(realpath "$(dirname "$0")/../..")
export PYTHONPATH=$REPO

python -m running_from_file.run_on_directory \
    --input-dir  "$REPO/Negatives/run_all_samples_through_newest_filter/data/filtered_full_old_wgermline" \
    --output-dir "$REPO/Negatives/run_TCGA_through_newest/results/with_fisher_fixed" \
    --workers 6
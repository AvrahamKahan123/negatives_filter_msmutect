#!/bin/bash
set -euo pipefail

#SRC_DIR=/home/avraham/MaruvkaLab/msmutect_postprocessing/run_msmutect_over_whole_tcga_again/
#DST_DIR=avrahamk@tech-ui01.hep.technion.ac.il:/storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/negatives_filter_msmutect/run_msmutect_over_whole_tcga_again

SRC_DIR=/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/complex_simulation/
DST_DIR=avrahamk@tech-ui01.hep.technion.ac.il:/storage/bfe_maruvka/avrahamk/positives_simulation/

# sync *.json, *.sh, *.py, and *.sub files from src_dir to dst_dir
rsync -avz --prune-empty-dirs \
    --include='*/' \
    --include='*.json' \
    --include='*.sh' \
    --include='*.py' \
    --include='*.sub' \
    --exclude='*' \
    "$@" \
    "$SRC_DIR" "$DST_DIR"

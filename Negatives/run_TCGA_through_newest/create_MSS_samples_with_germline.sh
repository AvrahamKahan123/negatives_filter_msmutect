#!/bin/bash

REPO=/home/avraham/MaruvkaLab/msmutect_postprocessing/
export PYTHONPATH=$REPO

python3 $REPO/converting_from_old_format_to_new/convert_old_msmutect_file_to_new.py --input-dir $REPO/data/mss_example_files \
      --output-dir data/mss_wgermline

python -m running_from_file.run_on_directory \
    --input-dir  data/mss_wgermline \
    --output-dir results/mss_samples_wfisher_fixed \
    --workers 6
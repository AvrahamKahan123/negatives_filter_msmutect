#!/bin/bash



export INPUT_FILE=/home/avraham/MaruvkaLab/msmutect_postprocessing/data/gib_files_filtered_oct1_2026/HG001_msmutect_normal0_tumor0.filt.mut.tsv
export ANNOTATED=annotated.filt.mut.tsv
export OUTPUT_FILE=tmp_tst

python -m converting_from_old_format_to_new.convert_old_msmutect_file_to_new \
        --input $INPUT_FILE --output $ANNOTATED
python3 -m running_from_file.run_single_sample --full-mut-tsv $ANNOTATED --output-prefix $OUTPUT_FILE
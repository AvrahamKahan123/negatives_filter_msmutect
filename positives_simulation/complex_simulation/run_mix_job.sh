#!/bin/bash
# Condor executable for the purity sweep: one sample, every purity.
#
#     run_mix_job.sh <sample>
#
# e.g. run_mix_job.sh TCGA-A6-5661
#
# mix_histograms.sub sets `getenv = true`, so the job inherits the environment you
# submitted from -- activate the conda environment with numpy BEFORE condor_submit and
# these jobs get it too.
#
# Every path is an environment variable with a default, so a trial run elsewhere needs no
# edit to this file:
#
#     OUTPUT_DIR=/storage/.../scratch condor_submit mix_histograms.sub

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <sample>" >&2
    exit 2
fi

SAMPLE=$1

# step 3 of run_msmutect_over_whole_tcga_again: current-MSMuTect calls, new format, and
# already cut down to CALL=="M" rows -- exactly the positives whose dilution we simulate
INPUT_DIR=${INPUT_DIR:-/storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/msmutect_results}
INPUT_SUFFIX=${INPUT_SUFFIX:-.full.mutated_only.mut.tsv}
OUTPUT_DIR=${OUTPUT_DIR:-/storage/bfe_maruvka/avrahamk/positives_simulation/mixed_histograms}
REPLICATES=${REPLICATES:-1}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON=${PYTHON:-python3}

INPUT="$INPUT_DIR/$SAMPLE$INPUT_SUFFIX"

echo "host=$(hostname) sample=$SAMPLE"
echo "python=$($PYTHON -c 'import sys; print(sys.executable)')"
echo "input=$INPUT"
echo "output=$OUTPUT_DIR/$SAMPLE"

if [ ! -e "$INPUT" ]; then
    echo "ERROR: $INPUT not found. Has step 3 of run_msmutect_over_whole_tcga_again run" >&2
    echo "       for this sample, and is INPUT_SUFFIX right?" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# FORCE=1 condor_submit ... redoes purities whose output already exists
ARGS=("$INPUT" "$OUTPUT_DIR/$SAMPLE" --replicates "$REPLICATES")
if [ "${FORCE:-0}" != "1" ]; then
    ARGS+=(--skip-existing)
fi
if [ -n "${SEED:-}" ]; then
    ARGS+=(--seed "$SEED")
fi

exec "$PYTHON" "$HERE/mix_histograms.py" "${ARGS[@]}"

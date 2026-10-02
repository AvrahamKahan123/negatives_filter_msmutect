#!/bin/bash
# Condor executable for all three steps.
#
#     run_job.sh <module> <sample>
#
# e.g. run_job.sh filter_to_mutations TCGA-W5-AA2K
#
# The submit files set `getenv = true`, so the job inherits the environment you submitted
# from -- activate MSMuTect's conda environment before condor_submit and these jobs get it
# too. PYTHONPATH is derived from this script's own location rather than hardcoded, so the
# checkout can live anywhere on the cluster.

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $0 <module> <sample>" >&2
    exit 2
fi

MODULE=$1
SAMPLE=$2

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(dirname "$HERE")
export PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}"

PYTHON=${PYTHON:-python3}

echo "host=$(hostname) sample=$SAMPLE module=$MODULE"
echo "python=$($PYTHON -c 'import sys; print(sys.executable)')"

# FORCE=1 condor_submit ... redoes samples whose output already exists
ARGS=(--sample "$SAMPLE")
if [ "${FORCE:-0}" = "1" ]; then
    ARGS+=(--force)
fi

exec "$PYTHON" -m "run_msmutect_over_whole_tcga_again.$MODULE" "${ARGS[@]}"

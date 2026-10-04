#!/bin/bash
# Submit the from_file sweep over every (sample, purity) pair.
#
#     conda activate genomics
#     ./submit_msmutect.sh
#
# Use this instead of calling `condor_submit msmutect_from_file.sub` directly. It
# rebuilds the job list from the histograms actually on disk, and creates the log
# directory -- condor does not create it, and holds every job instantly if it cannot open
# the log/output/error files.
#
# Any environment variable the submit file reads can be set here too, and is passed
# through by `getenv = true`:
#
#     RESULTS_DIR=/storage/.../scratch ./submit_msmutect.sh
#     FORCE=1 ./submit_msmutect.sh                 # redo pairs that already have output
#     SAMPLES="TCGA-A6-5661" ./submit_msmutect.sh  # just one sample

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# must match the default in msmutect_from_file.sub
LOG_DIR=${LOG_DIR:-/storage/bfe_maruvka/avrahamk/positives_simulation/logs/msmutect}
export LOG_DIR

echo "building job list..."
"$HERE/make_msmutect_job_list.sh"

mkdir -p "$LOG_DIR"
echo "logs -> $LOG_DIR"

# an executable that lost its +x in transit to the cluster is the other way every job
# gets held at once, and the hold reason does not say so in as many words
chmod +x "$HERE/run_msmutect_job.sh"

exec condor_submit "$HERE/msmutect_from_file.sub" "$@"

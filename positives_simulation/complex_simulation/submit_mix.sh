#!/bin/bash
# Submit the purity sweep.
#
#     ./submit_mix.sh
#
# Use this instead of calling `condor_submit mix_histograms.sub` directly.
#
# condor does NOT create the directories for log/output/error, and a job whose log files
# cannot be opened goes straight to Hold -- all of them, instantly, before anything runs.
# The sibling pipeline gets away with it because make_job_lists.py creates its log
# directory as a side effect of building the job list; this sweep queues from an inline
# list and so has no such step. Hence this wrapper.
#
# Any environment variable the submit file reads can be set here too, and is passed
# through by `getenv = true`:
#
#     OUTPUT_DIR=/storage/.../scratch ./submit_mix.sh
#     FORCE=1 ./submit_mix.sh

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# must match the default in mix_histograms.sub
LOG_DIR=${LOG_DIR:-/storage/bfe_maruvka/avrahamk/positives_simulation/logs/mix}
export LOG_DIR

mkdir -p "$LOG_DIR"
echo "logs -> $LOG_DIR"

# condor runs this itself, but an executable that lost its +x in transit is the other way
# every job gets held at once, and it is invisible in the hold reason until you look
chmod +x "$HERE/run_mix_job.sh"

exec condor_submit "$HERE/mix_histograms.sub" "$@"

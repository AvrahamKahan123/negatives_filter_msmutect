#!/bin/bash
# Condor executable for the from_file sweep: one (sample, purity) pair.
#
#     run_msmutect_job.sh <sample> <purity>
#
# e.g. run_msmutect_job.sh TCGA-A6-5661 0.15
#
# Calls mutations on the simulated pair written by the mixing step:
#   normal  <mixed_dir>/<sample>_1x.normal.hist.tsv          (one, shared by all purities)
#   tumor   <mixed_dir>/<sample>_1x_<purity>purity.hist.tsv
#
# msmutect.sh runs whichever python3 is on PATH, and the submit file sets getenv = true,
# so activate MSMuTect's conda environment BEFORE submitting.

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $0 <sample> <purity>" >&2
    exit 2
fi

SAMPLE=$1
PURITY=$2

MIXED_DIR=${MIXED_DIR:-/storage/bfe_maruvka/avrahamk/positives_simulation/mixed_histograms}
RESULTS_DIR=${RESULTS_DIR:-/storage/bfe_maruvka/avrahamk/positives_simulation/msmutect_results}
REPLICATES=${REPLICATES:-1}
MSMUTECT=${MSMUTECT:-/storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/MSMuTect_4/msmutect.sh}
READ_LEVEL=${READ_LEVEL:-5}

NORMAL="$MIXED_DIR/${SAMPLE}_${REPLICATES}x.normal.hist.tsv"
TUMOR="$MIXED_DIR/${SAMPLE}_${REPLICATES}x_${PURITY}purity.hist.tsv"
PREFIX="$RESULTS_DIR/${SAMPLE}_${REPLICATES}x_${PURITY}purity"
PRODUCED="$PREFIX.full.mut.tsv"

echo "host=$(hostname) sample=$SAMPLE purity=$PURITY"
echo "python=$(python3 -c 'import sys; print(sys.executable)')"
echo "normal=$NORMAL"
echo "tumor=$TUMOR"
echo "output=$PREFIX"

for f in "$NORMAL" "$TUMOR"; do
    if [ ! -e "$f" ]; then
        echo "ERROR: $f not found. Has the mixing step (./submit_mix.sh) finished?" >&2
        exit 1
    fi
done
if [ ! -x "$MSMUTECT" ]; then
    echo "ERROR: msmutect not executable at $MSMUTECT (set MSMUTECT=...)" >&2
    exit 1
fi

# FORCE=1 condor_submit ... redoes pairs that already have output
if [ -e "$PRODUCED" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "already called, skipping (FORCE=1 to redo): $PRODUCED"
    exit 0
fi

mkdir -p "$RESULTS_DIR"

# --from_file is single-core by design. -m because it is the only mode from_file supports;
# -f so a half-finished previous attempt at this pair is overwritten rather than refused.
"$MSMUTECT" --from_file -m -f -r "$READ_LEVEL" -N "$NORMAL" -T "$TUMOR" -O "$PREFIX"

# msmutect writes the tally itself; echo it so the count is in the job log too
COUNTS="$PREFIX.full.call_counts.tsv"
if [ -e "$COUNTS" ]; then
    echo "call counts:"
    cat "$COUNTS"
fi
echo "done -> $PRODUCED"

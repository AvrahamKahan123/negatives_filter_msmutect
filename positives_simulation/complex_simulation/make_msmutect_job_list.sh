#!/bin/bash
# Build the (sample, purity) job list for msmutect_from_file.sub.
#
#     ./make_msmutect_job_list.sh            # -> jobs/msmutect_from_file.txt
#
# The list is DISCOVERED from the mixed histograms on disk rather than regenerated from
# "0.1 to 1 in steps of 0.05". Two reasons:
#
#   1. The purity in the filename is whatever python's f"{round(x, 2)}" produced -- "0.1",
#      "0.95", "1.0". Re-deriving those strings in bash is a trap: printf "%g" turns 1.0
#      into "1", and every job for that purity would then look for a file that is not there.
#      Reading the names back cannot disagree with the names that exist.
#   2. A pair only makes it into the list if both of its files are actually present, so an
#      incomplete mixing run shows up as a short job list here instead of as a pile of
#      jobs that fail one by one on the cluster.
#
# Override the defaults with environment variables:
#   MIXED_DIR   where mix_histograms.sub wrote its output
#   REPLICATES  the Nx in the filenames (default 1)
#   SAMPLES     space-separated list; default is every sample found in MIXED_DIR

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

MIXED_DIR=${MIXED_DIR:-/storage/bfe_maruvka/avrahamk/positives_simulation/mixed_histograms}
REPLICATES=${REPLICATES:-1}
JOB_LIST=${JOB_LIST:-$HERE/jobs/msmutect_from_file.txt}

if [ ! -d "$MIXED_DIR" ]; then
    echo "ERROR: $MIXED_DIR does not exist. Has the mixing step (./submit_mix.sh) run?" >&2
    exit 1
fi

# a sample is one that has a normal histogram: that is the file all 19 purities pair with
if [ -z "${SAMPLES:-}" ]; then
    SAMPLES=""
    for normal in "$MIXED_DIR"/*_${REPLICATES}x.normal.hist.tsv; do
        [ -e "$normal" ] || continue
        name=$(basename "$normal")
        SAMPLES+=" ${name%_${REPLICATES}x.normal.hist.tsv}"
    done
fi

if [ -z "${SAMPLES// /}" ]; then
    echo "ERROR: no *_${REPLICATES}x.normal.hist.tsv in $MIXED_DIR." >&2
    echo "       Wrong REPLICATES, wrong MIXED_DIR, or the mixing step has not finished." >&2
    exit 1
fi

mkdir -p "$(dirname "$JOB_LIST")"
: > "$JOB_LIST"

total=0
for sample in $SAMPLES; do
    normal="$MIXED_DIR/${sample}_${REPLICATES}x.normal.hist.tsv"
    if [ ! -e "$normal" ]; then
        echo "  WARNING: $sample has no normal histogram, skipping all its purities" >&2
        continue
    fi
    count=0
    for tumor in "$MIXED_DIR/${sample}_${REPLICATES}x_"*purity.hist.tsv; do
        [ -e "$tumor" ] || continue
        purity=$(basename "$tumor")
        purity=${purity#${sample}_${REPLICATES}x_}
        purity=${purity%purity.hist.tsv}
        echo "$sample,$purity" >> "$JOB_LIST"
        count=$((count + 1))
    done
    printf "  %-16s %2d purities\n" "$sample" "$count"
    if [ "$count" -ne 19 ]; then
        echo "    NOTE: expected 19 (0.1 to 1.0 in steps of 0.05); mixing may be unfinished" >&2
    fi
    total=$((total + count))
done

echo "$total jobs -> $JOB_LIST"

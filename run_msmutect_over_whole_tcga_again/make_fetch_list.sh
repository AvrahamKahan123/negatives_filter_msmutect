#!/bin/bash
# Turn the unreadable-file report into the one-column list condor queues from.
#
#     ./make_fetch_list.sh [output_path]
#
# UNREAD_LIST may be "filename<TAB>reason" (what find_unreadable_files.py writes) or just
# filenames. Condor's `queue filename from <file>` assigns the WHOLE line to $(filename),
# so the reason column would end up glued onto the path. This takes the first field only.
#
# Files already recovered and verified in DEST_DIR are left out, so re-running after a
# partial batch gives exactly the work that remains.

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/fetch_config.sh"

OUT=${1:-$HERE/$FETCH_LIST}

if [ ! -r "$UNREAD_LIST" ]; then
    echo "ERROR: cannot read the unreadable-file list: $UNREAD_LIST" >&2
    echo "       Set UNREAD_LIST, or generate it with:" >&2
    echo "       python -m run_msmutect_over_whole_tcga_again.find_unreadable_files" >&2
    exit 1
fi

mkdir -p "$DEST_DIR"
# condor does NOT create its own log/output/error directories, and a missing one leaves the
# jobs held rather than running -- which looks exactly like "the jobs are hanging"
mkdir -p "$FETCH_LOG_DIR"

total=0
already=0
: > "$OUT"
while IFS= read -r line; do
    # first whitespace-separated field; skip blanks and comments
    name=$(printf '%s' "$line" | awk '{print $1}')
    [ -z "$name" ] && continue
    case "$name" in \#*) continue ;; esac
    name=$(basename "$name")          # tolerate a full path in the list
    total=$((total + 1))

    dest="$DEST_DIR/$name"
    if [ -f "$dest" ]; then
        # only count it as done if it is actually usable
        if [ "$VERIFY_GZIP" = "1" ] && [[ "$name" == *.gz ]]; then
            if gzip -t "$dest" 2>/dev/null; then already=$((already + 1)); continue; fi
            echo "  re-fetching (local copy fails gzip -t): $name" >&2
        else
            already=$((already + 1)); continue
        fi
    fi
    printf '%s\n' "$name" >> "$OUT"
done < "$UNREAD_LIST"

todo=$(wc -l < "$OUT" | tr -d ' ')
echo "from   $UNREAD_LIST"
echo "  $total file(s) listed"
echo "  $already already recovered in $DEST_DIR"
echo "  $todo to fetch -> $OUT"
[ "$todo" -eq 0 ] && echo "nothing to do."
exit 0

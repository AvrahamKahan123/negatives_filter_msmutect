#!/bin/bash
# Condor executable: recover ONE file from the backup host.
#
#     fetch_one_unread.sh TCGA-ZM-AA0N.full.mut.tsv.gz
#
# REMOTE_HOST:REMOTE_DIR/<filename>  ->  DEST_DIR/<filename>
#
# Downloads to a .partial and renames only after the copy succeeds AND (for .gz) the file
# decompresses, so an evicted or half-finished job can never leave something that looks
# recovered but is not. Re-running skips files already present and verified.

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/fetch_config.sh"

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <filename>" >&2
    exit 2
fi

NAME=$(basename "$1")
SRC="$REMOTE_DIR/$NAME"
DEST="$DEST_DIR/$NAME"
PARTIAL="$DEST.partial"

echo "host=$(hostname) file=$NAME"
echo "from=$REMOTE_HOST:$SRC"
echo "to=$DEST"

mkdir -p "$DEST_DIR"

verify() {
    # only .gz can be checked cheaply; anything else is accepted on a successful transfer
    if [ "$VERIFY_GZIP" = "1" ] && [[ "$NAME" == *.gz ]]; then
        gzip -t "$1"
    fi
}

# already done?
if [ -f "$DEST" ]; then
    if verify "$DEST" 2>/dev/null; then
        echo "already recovered and verified, skipping"
        exit 0
    fi
    echo "existing copy fails verification -- re-fetching"
fi

rm -f "$PARTIAL"

echo "scp opts: $SSH_OPTS ${SSH_KEY:+-i $SSH_KEY}"

# SSH_OPTS is deliberately unquoted so it word-splits into separate -o flags.
# BatchMode=yes matters: a condor job has no tty, so without it a password prompt would
# hang the job until eviction instead of failing.
# `|| rc=$?` rather than `if ! scp`: inside `if ! cmd` the then-branch sees $? == 0 (the
# status of the negation), so the real exit code is lost -- and a failed job would report
# success to condor.
rc=0
if [ "${FETCH_TIMEOUT:-0}" != "0" ] && command -v timeout >/dev/null 2>&1; then
    timeout --signal=TERM --kill-after=30 "$FETCH_TIMEOUT" \
        scp -p $SSH_OPTS ${SSH_KEY:+-i "$SSH_KEY"} "$REMOTE_HOST:$SRC" "$PARTIAL" || rc=$?
else
    scp -p $SSH_OPTS ${SSH_KEY:+-i "$SSH_KEY"} "$REMOTE_HOST:$SRC" "$PARTIAL" || rc=$?
fi
if [ "$rc" -eq 124 ]; then
    rm -f "$PARTIAL"
    echo "ERROR: $NAME timed out after ${FETCH_TIMEOUT}s -- the transfer stalled." >&2
    echo "       Raise FETCH_TIMEOUT if the file is genuinely huge or the link is slow." >&2
    exit 124
fi
if [ "$rc" -ne 0 ]; then
    rm -f "$PARTIAL"
    echo "ERROR: scp failed (exit $rc) for $NAME" >&2
    echo "       If this was a permission/publickey error: condor jobs cannot type a" >&2
    echo "       password. Set up key-based auth to $REMOTE_HOST and point SSH_KEY at the" >&2
    echo "       key, or run the fetch from the login node instead of condor." >&2
    exit "$rc"
fi

if ! verify "$PARTIAL"; then
    rm -f "$PARTIAL"
    echo "ERROR: $NAME transferred but does not decompress -- the backup copy is damaged too." >&2
    exit 1
fi

mv -f "$PARTIAL" "$DEST"
echo "recovered $(stat -c%s "$DEST" 2>/dev/null || echo '?') bytes -> $DEST"

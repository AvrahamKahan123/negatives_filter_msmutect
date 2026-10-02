#!/bin/bash
# Keep only called mutations (CALL == "M") from a .full.mut.tsv.gz file.
#
#   usage: ./filter_called_mutations.sh <input.full.mut.tsv.gz> [CALL_COLUMN]
#
# Writes <input>.filt.mut.tsv.gz next to the input, header preserved.
#
# The CALL column is located by name in the header, because its index is not
# stable across msmutect versions: it is column 51 in the 63 field schema
# (Texas/entropy_analysis) but column 50 in the 57 field schema
# (msmutect_development). Pass CALL_COLUMN to override for headerless files.
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "usage: $0 <input.full.mut.tsv.gz> [CALL_COLUMN]" >&2
    exit 1
fi

IN="$1"
FORCED_COL="${2:-0}"

case "$IN" in
    *.full.mut.tsv.gz) ;;
    *) echo "error: expected a .full.mut.tsv.gz file, got: $IN" >&2; exit 1 ;;
esac

if [ ! -r "$IN" ]; then
    echo "error: cannot read $IN" >&2
    exit 1
fi

OUT="${IN%.full.mut.tsv.gz}.filt.mut.tsv.gz"

if [ -e "$OUT" ]; then
    echo "error: $OUT already exists, refusing to overwrite" >&2
    exit 1
fi

# Write to a temp file first so a failure never leaves a truncated .gz behind.
TMP="$(mktemp "${OUT}.partial.XXXXXX")"
trap 'rm -f "$TMP"' EXIT

# Single pass: resolve the CALL column from the header, then filter.
# Match is exact, so "NM" is not mistaken for "M".
zcat -- "$IN" \
  | awk -F'\t' -v forced="$FORCED_COL" '
      NR == 1 {
          if (forced > 0) {
              call = forced
          } else {
              for (i = 1; i <= NF; i++)
                  if ($i == "CALL") call = i
              if (!call) {
                  print "error: no CALL column in header; pass the index explicitly" > "/dev/stderr"
                  exit 1
              }
          }
          print
          next
      }
      $call == "M" { kept++; print }
      END { printf "kept %d mutation lines (CALL column %d)\n", kept, call > "/dev/stderr" }
  ' \
  | gzip > "$TMP"

mv "$TMP" "$OUT"
trap - EXIT

echo "wrote $OUT"

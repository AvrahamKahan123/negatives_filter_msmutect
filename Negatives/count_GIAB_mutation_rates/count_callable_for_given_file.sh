#!/bin/bash

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <remote_path>" >&2
    exit 1
fi

REMOTE_PATH="$1"
LOCAL_FILE="$(basename "$REMOTE_PATH")"

scp "avrahamk@tech-bfe-st01:${REMOTE_PATH}" "$LOCAL_FILE"
python3 count_callable.py "$LOCAL_FILE"
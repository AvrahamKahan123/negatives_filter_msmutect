#!/bin/bash
# Shared settings for the unreadable-file recovery (make_fetch_list.sh / fetch_one_unread.sh).
# Every value can be overridden from the environment before condor_submit, e.g.
#
#     FETCH_WORKERS=4 condor_submit fetch_unread.sub
#
# Nothing here is specific to a single file; the per-file work is in fetch_one_unread.sh.

# where the good copies live (the BACKUP tree, not the one throwing I/O errors)
: "${REMOTE_HOST:=avrahamk@tech-bfe-st01}"
: "${REMOTE_DIR:=/storage/bfe_maruvka/gaiafr/Folders_backup/Research_project/WGS_SNVs_indels_analysis_project/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run}"

# where recovered files land
: "${DEST_DIR:=/storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/unread}"

# the list produced by find_unreadable_files.py (or written by hand)
: "${UNREAD_LIST:=/storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/negatives_filter_msmutect/run_msmutect_over_whole_tcga_again/jobs/unread_files.txt}"

# one filename per line, derived from UNREAD_LIST by make_fetch_list.sh. Condor's
# `queue ... from` hands the WHOLE line to the job, so the list it reads must be one clean
# column -- UNREAD_LIST may carry a tab-separated reason alongside each name.
: "${FETCH_LIST:=fetch_list.txt}"

# ssh must be non-interactive: a condor job has no tty, so a password prompt would hang the
# job until it is evicted rather than fail. BatchMode makes it exit immediately instead.
: "${SSH_KEY:=}"                      # optional: path to the private key to use
: "${SSH_OPTS:=-o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=30}"

# verify each recovered .gz actually decompresses. These are replacements for files the
# filesystem could not read, so "it transferred" is not the same as "it is good".
: "${VERIFY_GZIP:=1}"

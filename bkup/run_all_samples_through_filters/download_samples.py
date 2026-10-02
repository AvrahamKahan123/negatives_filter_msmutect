import argparse
import os

import config
import ssh_connection

# Step 1 (see Plan.md): scp/rsync the 5 samples' filtered and full MSMuTect output files
# from the lab server to local subdirectories of this directory.
#
# Auth note: no SSH key from this machine works against REMOTE_HOST, so you have to type
# your password -- but only once. Every transfer reuses a single shared SSH connection
# (see ssh_connection.py), so run this from a real terminal where the prompt can reach you.


def rsync(remote_path: str, local_path: str, force: bool):
    if os.path.exists(local_path) and not force:
        print(f"already have {os.path.basename(local_path)}, skipping (use --force to re-download)")
        return
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    print(f"downloading {remote_path} -> {local_path}")
    result = ssh_connection.rsync_over_master(remote_path, local_path)
    if result.returncode != 0:
        raise RuntimeError(f"rsync failed ({result.returncode}) for {remote_path}")


def download_filt(sample: str, force: bool = False):
    config.ensure_free_disk_gb()
    rsync(config.remote_filt_path(sample), config.filt_gz_local_path(sample), force)


def download_full(sample: str, force: bool = False):
    config.ensure_free_disk_gb()
    rsync(config.remote_full_path(sample), config.full_gz_local_path(sample), force)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES,
                         help="subset of samples to download (default: all 5 in config.SAMPLES)")
    parser.add_argument("--full", action="store_true",
                         help="also download the (large) *.full.mut.tsv.gz files for every sample "
                              "up front, instead of letting rerun_msmutect_from_file.py fetch them "
                              "one at a time as it processes each sample")
    parser.add_argument("--force", action="store_true", help="re-download even if the local file already exists")
    args = parser.parse_args()

    for sample in args.samples:
        download_filt(sample, force=args.force)
        if args.full:
            download_full(sample, force=args.force)


if __name__ == "__main__":
    main()

import argparse
import os

import config
import download_samples
import ssh_connection

# Fetches EVERY input the pipeline needs (both *.called.filt.mut.tsv.gz and
# *.full.mut.tsv.gz for all 5 samples) before any analysis starts, then hands the SSH
# connection back.
#
# Why up front: a single sample's MSMuTect --from_file run takes hours, so downloading the
# next sample's file lazily in between would come back to a control master that has since
# expired (ControlPersist=4h) -- either re-prompting for a password or, in a nohup'd run
# with no terminal, failing outright. Downloading everything first means the whole rest of
# the pipeline is offline, and run_parallel.py can start 5 samples at once without any of
# them racing for the connection.
#
# Expect roughly 0.5 GB per *.full.mut.tsv.gz plus ~30-110 MB per filt file.


def total_bytes(paths):
    return sum(os.path.getsize(p) for p in paths if os.path.exists(p))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--force", action="store_true", help="re-download even if a local copy exists")
    parser.add_argument("--keep-ssh", action="store_true",
                         help="leave the shared SSH connection open when finished "
                              "(by default it is closed, since nothing later needs it)")
    args = parser.parse_args()

    for i, sample in enumerate(args.samples, 1):
        print(f"\n===== [{i}/{len(args.samples)}] {sample} =====")
        download_samples.download_filt(sample, force=args.force)
        download_samples.download_full(sample, force=args.force)

    filt_paths = [config.filt_gz_local_path(s) for s in args.samples]
    full_paths = [config.full_gz_local_path(s) for s in args.samples]
    missing = [p for p in filt_paths + full_paths if not os.path.exists(p)]

    print("\n===== downloaded =====")
    print(f"filt files: {len(filt_paths) - len([p for p in filt_paths if not os.path.exists(p)])}"
          f"/{len(filt_paths)}  ({total_bytes(filt_paths) / 1024 ** 3:.2f} GB)")
    print(f"full files: {len(full_paths) - len([p for p in full_paths if not os.path.exists(p)])}"
          f"/{len(full_paths)}  ({total_bytes(full_paths) / 1024 ** 3:.2f} GB)")

    if not args.keep_ssh:
        ssh_connection.close_master()

    if missing:
        raise SystemExit("MISSING after download:\n  " + "\n  ".join(missing))
    print("\nall inputs present -- the rest of the pipeline needs no network access.")
    print("next: python run_parallel.py")


if __name__ == "__main__":
    main()

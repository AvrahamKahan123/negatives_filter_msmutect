import argparse
import os
import subprocess
import sys
import time

import config
import compare_old_vs_new

# Runs the whole pipeline on every sample AT THE SAME TIME, one process per sample.
#
# MSMuTect's --from_file mode is single-core by design, so a sample cannot be sped up
# internally -- the way to use this machine's cores is to give each sample its own process.
# Each worker is just
#     python run_all.py --samples <SAMPLE> --skip-summary
# so all the actual work is the existing per-sample code path; this file only supervises.
#
# Every worker writes only its own sample's files, and step 4 saves a per-sample row rather
# than a shared summary, so the workers never write to the same path. summary.tsv is built
# here once they have all finished.
#
# Run download_all.py FIRST. If the inputs are already local the workers never touch the
# network, so five of them starting at once cannot contend for (or outlive) the SSH
# connection. This script refuses to start if anything is still missing.


def missing_inputs(samples):
    missing = []
    for sample in samples:
        for path in (config.filt_gz_local_path(sample), config.full_gz_local_path(sample)):
            if not os.path.exists(path):
                missing.append(path)
    return missing


def launch(sample: str, force: bool) -> tuple:
    os.makedirs(config.LOGS_DIR, exist_ok=True)
    log_path = config.sample_log_path(sample)
    command = [sys.executable, os.path.join(config.THIS_DIR, "run_all.py"),
               "--samples", sample, "--skip-summary"]
    if force:
        command.append("--force")
    log_file = open(log_path, "w")
    # line-buffered through the log so `tail -f` shows progress as it happens
    process = subprocess.Popen(command, cwd=config.THIS_DIR, stdout=log_file,
                               stderr=subprocess.STDOUT, text=True)
    print(f"  {sample:<16} pid {process.pid:<8} -> {log_path}")
    return sample, process, log_file


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--force", action="store_true", help="redo every step even if outputs already exist")
    parser.add_argument("--max-workers", type=int, default=None,
                         help="how many samples to run at once (default: all of them, since "
                              "each is a single core and there are only 5)")
    parser.add_argument("--allow-missing-inputs", action="store_true",
                         help="start even if some inputs are not downloaded yet; workers will "
                              "then fetch them over SSH themselves")
    args = parser.parse_args()

    if not args.allow_missing_inputs:
        missing = missing_inputs(args.samples)
        if missing:
            raise SystemExit(
                f"{len(missing)} input file(s) not downloaded yet:\n  " + "\n  ".join(missing) +
                "\n\nRun `python download_all.py` first so the workers need no network access "
                "(or pass --allow-missing-inputs to let them download as they go).")

    max_workers = args.max_workers or len(args.samples)
    config.ensure_free_disk_gb()
    print(f"starting {len(args.samples)} sample(s), {max_workers} at a time")

    started = time.time()
    pending = list(args.samples)
    running = []
    results = {}

    while pending or running:
        while pending and len(running) < max_workers:
            running.append(launch(pending.pop(0), args.force))

        time.sleep(5)
        for entry in list(running):
            sample, process, log_file = entry
            if process.poll() is None:
                continue
            running.remove(entry)
            log_file.close()
            results[sample] = process.returncode
            elapsed = (time.time() - started) / 3600
            status = "ok" if process.returncode == 0 else f"FAILED (exit {process.returncode})"
            print(f"  {sample:<16} {status}  [{elapsed:.2f} h elapsed, "
                  f"{len(running)} running, {len(pending)} queued]")

    print(f"\nall workers finished in {(time.time() - started) / 3600:.2f} h")
    failed = [s for s, code in results.items() if code != 0]
    for sample in failed:
        print(f"  {sample} failed -- see {config.sample_log_path(sample)}")

    print("\n===== summary =====")
    compare_old_vs_new.summarize(args.samples)

    if failed:
        raise SystemExit(f"{len(failed)} sample(s) failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()

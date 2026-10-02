"""Find every file in a directory that the filesystem cannot serve.

Scans a directory and records which files fail, and HOW -- separating the storage-level
failures (Input/output error, stale handle, permissions) from files that are merely absent.
`ls` reports only the first one it hits and stops being useful across thousands of files;
this walks the whole directory and gives you a list plus a breakdown by errno.

    # the whole-TCGA source directory, as `ls` sees it
    python -m run_msmutect_over_whole_tcga_again.find_unreadable_files

    # any directory, and also try actually reading each file's first bytes
    python -m run_msmutect_over_whole_tcga_again.find_unreadable_files --dir /some/path --read

    # is the set stable, or is the mount flaky?
    python -m run_msmutect_over_whole_tcga_again.find_unreadable_files --repeat 3

    # if /storage is Lustre: which OST do the bad files live on?
    python -m run_msmutect_over_whole_tcga_again.find_unreadable_files --lfs

Checking levels, cheapest first:
    (default)  os.stat  -- exactly what `ls` does; catches metadata-level faults
    --read     stat + read the first block -- catches a disk that fails on data but not inodes
    --gzip     full decompression -- also catches genuinely corrupt CONTENT, but reads every
               byte of every file, so only use it on a short list

Note that corrupt content and an I/O error are different layers: a truncated or garbled .gz
stats perfectly and only fails when decompressed, so --gzip is the only level that sees it.
"""

import argparse
import errno
import fnmatch
import gzip
import os
import subprocess
import sys
from collections import Counter, defaultdict

from run_msmutect_over_whole_tcga_again import config

READ_BLOCK = 1 << 20  # 1 MiB
OK = "ok"


def classify(path: str, check_read: bool = False, check_gzip: bool = False) -> str:
    """-> OK, or a short reason string naming the errno."""
    try:
        os.stat(path)
    except OSError as e:
        return f"stat {errno.errorcode.get(e.errno, e.errno)}: {e.strerror}"

    if check_gzip:
        try:
            with gzip.open(path, "rb") as f:
                while f.read(READ_BLOCK):
                    pass
        except OSError as e:
            code = errno.errorcode.get(getattr(e, "errno", None), "")
            return f"read {code}: {e.strerror or e}".strip()
        except Exception as e:  # zlib/BadGzipFile: content is damaged, not the storage
            return f"corrupt content: {type(e).__name__}: {e}"
        return OK

    if check_read:
        try:
            with open(path, "rb") as f:
                f.read(READ_BLOCK)
        except OSError as e:
            return f"read {errno.errorcode.get(e.errno, e.errno)}: {e.strerror}"
    return OK


def list_directory(directory: str, pattern: str) -> list:
    """Names only -- os.listdir does not stat, so one bad file cannot break the listing."""
    try:
        names = os.listdir(directory)
    except OSError as e:
        raise SystemExit(f"ERROR: cannot list {directory}: {e.strerror}\n"
                         f"       The directory itself is unreadable, not just files in it.")
    return sorted(n for n in names if fnmatch.fnmatch(n, pattern))


def lustre_ost(path: str) -> str:
    """Which Lustre OST holds this file, if lfs is available. Layout lives on the MDT, so
    this often works even when the file's DATA cannot be read -- which is the whole point."""
    try:
        out = subprocess.run(["lfs", "getstripe", path], capture_output=True, text=True,
                             timeout=30)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return ""
    if out.returncode != 0:
        return ""
    indices = []
    seen_header = False
    for line in out.stdout.splitlines():
        if "obdidx" in line:
            seen_header = True
            continue
        if seen_header and line.strip():
            first = line.split()[0]
            if first.isdigit():
                indices.append(first)
    return ",".join(indices)


def scan(directory: str, pattern: str, check_read: bool, check_gzip: bool,
         show_progress: bool = True) -> dict:
    names = list_directory(directory, pattern)
    if not names:
        raise SystemExit(f"no files matching {pattern!r} in {directory}")
    results = {}
    for i, name in enumerate(names, 1):
        results[name] = classify(os.path.join(directory, name), check_read, check_gzip)
        if show_progress and i % 500 == 0:
            bad = sum(1 for r in results.values() if r != OK)
            print(f"  checked {i:,}/{len(names):,} ({bad:,} bad so far)", flush=True)
    return results


def report(results: dict, directory: str, out_path: str, use_lfs: bool):
    bad = {name: reason for name, reason in results.items() if reason != OK}
    print(f"\n{len(results):,} file(s) checked in {directory}")
    print(f"  readable   : {len(results) - len(bad):,}")
    print(f"  UNREADABLE : {len(bad):,}")
    if not bad:
        print("\nall files are readable.")
        return bad

    print("\nby failure:")
    for reason, count in Counter(bad.values()).most_common():
        print(f"  {count:>7,}  {reason}")

    print("\nfirst few:")
    for name in sorted(bad)[:10]:
        print(f"  {name}  ({bad[name]})")
    if len(bad) > 10:
        print(f"  ... and {len(bad) - 10:,} more")

    osts = defaultdict(list)
    if use_lfs:
        print("\nlooking up Lustre layout for the bad files...")
        for name in sorted(bad):
            osts[lustre_ost(os.path.join(directory, name)) or "unknown"].append(name)
        if list(osts) != ["unknown"]:
            print("bad files per OST:")
            for ost, names in sorted(osts.items(), key=lambda kv: -len(kv[1])):
                print(f"  OST {ost:<12} {len(names):>7,} file(s)")
            print("\n  If they concentrate on one OST, that is almost certainly a single")
            print("  failing storage target -- give the admin that OST number.")
        else:
            print("  (lfs returned nothing useful; this may not be Lustre)")

    with open(out_path, "w") as f:
        for name in sorted(bad):
            line = f"{name}\t{bad[name]}"
            if use_lfs:
                ost = next((o for o, names in osts.items() if name in names), "")
                line += f"\tOST={ost}"
            f.write(line + "\n")
    print(f"\nfull list -> {out_path}")
    return bad


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", default=config.INPUT_DIR, help="directory to scan")
    parser.add_argument("--pattern", default="*" + config.INPUT_SUFFIX, help="filename glob")
    parser.add_argument("--read", action="store_true",
                        help="also read the first 1 MiB of each file")
    parser.add_argument("--gzip", action="store_true",
                        help="fully decompress each file (slow; also finds corrupt CONTENT)")
    parser.add_argument("--lfs", action="store_true",
                        help="report which Lustre OST each bad file is on")
    parser.add_argument("--repeat", type=int, default=1,
                        help="scan N times and report whether the bad set is stable")
    parser.add_argument("--output", default=None, help="where to write the bad-file list")
    args = parser.parse_args()

    out_path = args.output or os.path.join(config.JOBS_DIR, "unreadable_files.txt")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    print(f"scanning {args.dir}")
    print(f"  pattern {args.pattern}")
    print(f"  level   {'full gzip' if args.gzip else 'stat + read' if args.read else 'stat (same as ls)'}")

    runs = []
    for attempt in range(1, args.repeat + 1):
        if args.repeat > 1:
            print(f"\n--- pass {attempt}/{args.repeat} ---")
        results = scan(args.dir, args.pattern, args.read, args.gzip)
        runs.append({n for n, r in results.items() if r != OK})

    bad = report(results, args.dir, out_path, args.lfs)

    if args.repeat > 1:
        always = set.intersection(*runs)
        ever = set.union(*runs)
        print(f"\nacross {args.repeat} passes: {len(always):,} failed every time, "
              f"{len(ever):,} failed at least once")
        if len(always) != len(ever):
            print("  The set CHANGES between passes -- that points at a flaky mount or an")
            print("  intermittent storage fault rather than specific damaged files.")
        else:
            print("  The same files fail every time -- consistent with specific bad files or")
            print("  one consistently unreachable storage target.")

    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

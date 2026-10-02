import argparse
import csv
import os
import subprocess
from collections import Counter

import pandas as pd

import config
import old_mutation_calls
from Negatives.run_all_samples_through_newest_filter import rerun_msmutect_from_file as rerun

# Measures what switching call_verified_locus to reversion_to_reference_strict actually does,
# WITHOUT re-running MSMuTect over all 27.7M loci.
#
# Every locus fed in here was called M (Mutation) by the previous --from_file run, so any
# locus that now comes back as something else (RR, in practice) is a locus the change moved.
# Two groups are run together so the win and the cost are measured on the same run:
#
#   violations : new-version mutations ABSENT from the old version's call set
#                (results/yossi_encapsulation/*.new_not_in_old_raw.tsv). Every one of these
#                that flips to RR is a subset violation repaired.
#   control    : an evenly-spaced sample of the new version's OTHER mutations, which ARE in
#                the old call set. Every one of these that flips to RR is collateral loss.
#
# Histograms are carved by streaming the full file and the loci file in lockstep -- same
# alignment check as rerun_msmutect_from_file -- and emitting only the selected loci, so the
# germline annotation stays correct even though the output is a sparse subset.

EXPERIMENT_DIR = os.path.join(config.RESULTS_DIR, "strict_rr_experiment")
VIOLATION, CONTROL = "violation", "control"


def violation_loci(sample: str) -> set:
    fp = os.path.join(config.ENCAPSULATION_DIR, f"{sample}.new_not_in_old_raw.tsv")
    if not os.path.exists(fp):
        raise FileNotFoundError(f"{fp} not found -- the new-vs-old raw comparison has not been run")
    df = pd.read_csv(fp, delimiter="\t", low_memory=False)
    return {tuple(str(v) for v in r) for r in df[config.LOCUS_KEY_COLUMNS].itertuples(index=False, name=None)}


def control_loci(sample: str, violations: set, size: int) -> set:
    # evenly spaced rather than random so the selection is reproducible
    df = pd.read_csv(config.new_mutated_only_path(sample), delimiter="\t",
                     low_memory=False, usecols=config.LOCUS_KEY_COLUMNS)
    everything = [tuple(str(v) for v in r)
                  for r in df[config.LOCUS_KEY_COLUMNS].itertuples(index=False, name=None)]
    eligible = [k for k in everything if k not in violations]
    if size >= len(eligible):
        return set(eligible)
    stride = len(eligible) / size
    return {eligible[int(i * stride)] for i in range(size)}


def carve_subset(sample: str, wanted: set, normal_dest: str, tumor_dest: str) -> int:
    """Carve histograms for `wanted` only, taking germline line-for-line from the loci file."""
    source = config.full_gz_local_path(sample)
    if not os.path.exists(source):
        source = source[:-3]  # the .gz may have been decompressed in place
    if not os.path.exists(source):
        raise FileNotFoundError(f"no full file for {sample} in {config.FULL_DIR}")

    written = 0
    with old_mutation_calls.open_maybe_gzip(source) as src, \
            open(config.GERMLINE_LOCI_FILE) as loci, \
            open(normal_dest, "w", newline="\n") as n_out, \
            open(tumor_dest, "w", newline="\n") as t_out:
        reader = csv.reader(src, delimiter="\t")
        n_writer = csv.writer(n_out, delimiter="\t", lineterminator="\n")
        t_writer = csv.writer(t_out, delimiter="\t", lineterminator="\n")
        header = next(reader)
        index = {name: i for i, name in enumerate(header)}
        n_plan = rerun.column_plan(index, config.LOCUS_COLUMNS + rerun.histogram_columns("NORMAL_"))
        t_plan = rerun.column_plan(index, config.LOCUS_COLUMNS + rerun.histogram_columns("TUMOR_"))
        n_writer.writerow([n for n, _ in n_plan])
        t_writer.writerow([n for n, _ in t_plan])
        key_idx = [index[c] for c in config.LOCUS_KEY_COLUMNS]

        for row_number, row in enumerate(reader, 1):
            loci_line = loci.readline()
            if not loci_line:
                raise RuntimeError(f"loci file ran out at row {row_number} of {source}")
            loci_locus, annotation = rerun.parse_loci_line(loci_line)
            mut_locus = tuple(row[i] for i in key_idx)
            if mut_locus[:3] != loci_locus[:3]:
                raise RuntimeError(f"locus mismatch at row {row_number}: {mut_locus} vs {loci_locus}")
            if mut_locus not in wanted:
                continue
            germline = "True" if annotation == rerun.HETEROGENEOUS else "False"
            n_writer.writerow([germline if i is None else row[i] for _, i in n_plan])
            t_writer.writerow([germline if i is None else row[i] for _, i in t_plan])
            written += 1
            if written == len(wanted):
                break
    return written


def run_msmutect(normal_hist: str, tumor_hist: str, output_prefix: str) -> str:
    command = [config.MSMUTECT_PATH, "-N", normal_hist, "-T", tumor_hist,
               "-m", "-A", "-O", output_prefix, "--from_file", "-f"]
    print("  " + " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"msmutect failed ({result.returncode})\n"
                           f"stdout:\n{result.stdout[-3000:]}\nstderr:\n{result.stderr[-3000:]}")
    return output_prefix + ".full.mut.tsv"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", default="TCGA-A6-5661")
    parser.add_argument("--control-size", type=int, default=20000)
    parser.add_argument("--reuse-histograms", action="store_true",
                         help="skip carving if the histogram files are already there")
    args = parser.parse_args()
    sample = args.sample

    os.makedirs(EXPERIMENT_DIR, exist_ok=True)
    violations = violation_loci(sample)
    controls = control_loci(sample, violations, args.control_size)
    wanted = violations | controls
    print(f"[{sample}] {len(violations):,} violation loci + {len(controls):,} control loci "
          f"= {len(wanted):,} loci (all were M in the previous run)")

    normal_hist = os.path.join(EXPERIMENT_DIR, f"{sample}.subset_normal.hist.tsv")
    tumor_hist = os.path.join(EXPERIMENT_DIR, f"{sample}.subset_tumor.hist.tsv")
    if args.reuse_histograms and os.path.exists(normal_hist) and os.path.exists(tumor_hist):
        print(f"[{sample}] reusing existing histograms")
    else:
        print(f"[{sample}] carving subset histograms (streams the full file + loci file)")
        found = carve_subset(sample, wanted, normal_hist, tumor_hist)
        print(f"[{sample}] carved {found:,}/{len(wanted):,} loci")

    print(f"[{sample}] running msmutect --from_file with the CURRENT CallMutations.py")
    output = run_msmutect(normal_hist, tumor_hist, os.path.join(EXPERIMENT_DIR, f"{sample}.subset"))

    calls = {}
    with open(output, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            calls[tuple(row[c] for c in config.LOCUS_KEY_COLUMNS)] = row["CALL"]

    rows = []
    for label, group in ((VIOLATION, violations), (CONTROL, controls)):
        counts = Counter(calls.get(k, "<missing>") for k in group)
        still_m = counts.get("M", 0)
        moved = len(group) - still_m
        rows.append({"GROUP": label, "LOCI": len(group), "STILL_M": still_m,
                     "MOVED_OFF_M": moved,
                     "MOVED_PCT": round(100 * moved / max(len(group), 1), 2),
                     **{f"NOW_{k}": v for k, v in counts.most_common()}})

    summary = pd.DataFrame(rows)
    summary_fp = os.path.join(EXPERIMENT_DIR, f"{sample}.strict_rr_effect.tsv")
    summary.to_csv(summary_fp, sep="\t", index=False)
    print()
    print(summary.to_string(index=False))
    print(f"\n-> {summary_fp}")

    v = rows[0]
    c = rows[1]
    print(f"\n[{sample}] subset violations repaired : {v['MOVED_OFF_M']:,}/{v['LOCI']:,} "
          f"({v['MOVED_PCT']}%)")
    print(f"[{sample}] control mutations lost      : {c['MOVED_OFF_M']:,}/{c['LOCI']:,} "
          f"({c['MOVED_PCT']}%)")
    total_new = len(pd.read_csv(config.new_mutated_only_path(sample), delimiter="\t",
                                low_memory=False, usecols=["CHROMOSOME"]))
    non_viol = total_new - len(violations)
    print(f"[{sample}] extrapolated loss over all {non_viol:,} in-old mutations: "
          f"~{round(non_viol * c['MOVED_PCT'] / 100):,}")


if __name__ == "__main__":
    main()

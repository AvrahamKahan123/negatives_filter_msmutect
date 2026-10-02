import argparse
import json
import os

import pandas as pd

import config
import msmutect_format
import old_mutation_calls
from results_postprocessing.analyze_full_tcga_mut_set import yossi_filter

# Step 4: verify that the newest MSMuTect marks as mutations a SUBSET of what the old
# version did (this is the actual goal stated at the top of Plan.md). Two comparisons are
# made per sample, both keyed on locus = (CHROMOSOME, START, END, PATTERN):
#
#   raw:   new MSMuTect's own calls (--from_file rerun's *.full.mutated_only.mut.tsv)
#          vs. the old run's own CALL=="M" rows (old_mutation_calls.py).
#   yossi: the current yossi_filter() applied to each of the above.
#
# Both sides of the raw comparison are each version's UNFILTERED call set, taken from the
# respective full runs. *.called.filt.mut.tsv.gz is deliberately not used as the baseline:
# it is already filtered and omits some of the old version's calls, which would show up
# here as new-version violations.
#
# New-format frames go through msmutect_format.normalize_for_yossi_filter first -- current
# MSMuTect renamed FISHER_TEST_P_VALUE and now writes empty (not 0) allele slots, either of
# which makes yossi_filter fail or silently return almost nothing.
#
# Requires run_yossi_filter_on_old.py and rerun_msmutect_from_file.py to have already run.


def locus_keys(df: pd.DataFrame) -> set:
    return set(df[config.LOCUS_KEY_COLUMNS].itertuples(index=False, name=None))


def read_call_counts(sample: str) -> dict:
    fp = config.rerun_output_prefix(sample) + ".full.call_counts.tsv"
    if not os.path.exists(fp):
        return {}
    counts_df = pd.read_csv(fp, sep="\t")
    return dict(zip(counts_df["CALL"], counts_df["count"]))


def new_yossi_filtered(sample: str, noisy_db) -> pd.DataFrame:
    mutated_only_fp = config.rerun_output_prefix(sample) + ".full.mutated_only.mut.tsv"
    if not os.path.exists(mutated_only_fp):
        raise FileNotFoundError(f"{mutated_only_fp} not found -- run rerun_msmutect_from_file.py first")
    new_mutations_df = pd.read_csv(mutated_only_fp, delimiter="\t", low_memory=False)
    normalized = msmutect_format.normalize_for_yossi_filter(new_mutations_df)
    print(f"[{sample}] new-format normalization: {msmutect_format.describe_normalization(new_mutations_df)}")
    filtered = yossi_filter(normalized, noisy_db, name=sample)
    os.makedirs(config.NEW_YOSSI_DIR, exist_ok=True)
    filtered.to_csv(config.new_yossi_filtered_path(sample), sep="\t", index=False)
    return new_mutations_df, filtered


def compare_one(sample: str, noisy_db=None) -> dict:
    print(f"[{sample}] comparing old vs. new")

    old_mutations_df = pd.read_csv(old_mutation_calls.extract(sample), delimiter="\t", low_memory=False)
    old_yossi_fp = config.old_yossi_filtered_path(sample)
    if not os.path.exists(old_yossi_fp):
        raise FileNotFoundError(f"{old_yossi_fp} not found -- run run_yossi_filter_on_old.py first")
    old_yossi_df = pd.read_csv(old_yossi_fp, delimiter="\t", low_memory=False)

    new_mutations_df, new_yossi_df = new_yossi_filtered(sample, noisy_db or config.noisy_locus_db())

    old_mut_keys = locus_keys(old_mutations_df)
    new_mut_keys = locus_keys(new_mutations_df)
    old_yossi_keys = locus_keys(old_yossi_df)
    new_yossi_keys = locus_keys(new_yossi_df)

    raw_violations = new_mut_keys - old_mut_keys
    yossi_violations = new_yossi_keys - old_yossi_keys

    os.makedirs(config.COMPARISON_DIR, exist_ok=True)
    if raw_violations:
        pd.DataFrame(sorted(raw_violations), columns=config.LOCUS_KEY_COLUMNS).to_csv(
            os.path.join(config.COMPARISON_DIR, f"{sample}.raw_violations.tsv"), sep="\t", index=False)
    if yossi_violations:
        pd.DataFrame(sorted(yossi_violations), columns=config.LOCUS_KEY_COLUMNS).to_csv(
            os.path.join(config.COMPARISON_DIR, f"{sample}.yossi_violations.tsv"), sep="\t", index=False)

    call_counts = read_call_counts(sample)
    row = {
        "SAMPLE": sample,
        "OLD_MUTATIONS": len(old_mut_keys),
        "NEW_MUTATIONS": len(new_mut_keys),
        "RAW_SUBSET_OK": len(raw_violations) == 0,
        "RAW_VIOLATIONS": len(raw_violations),
        "OLD_YOSSI_FILTERED": len(old_yossi_keys),
        "NEW_YOSSI_FILTERED": len(new_yossi_keys),
        "YOSSI_SUBSET_OK": len(yossi_violations) == 0,
        "YOSSI_VIOLATIONS": len(yossi_violations),
    }
    row.update({f"NEW_CALL_{k}": v for k, v in call_counts.items()})

    status = "PASS" if row["RAW_SUBSET_OK"] and row["YOSSI_SUBSET_OK"] else "FAIL"
    print(f"[{sample}] {status}: old={row['OLD_MUTATIONS']} new={row['NEW_MUTATIONS']} "
          f"raw_violations={row['RAW_VIOLATIONS']} | "
          f"old_yossi={row['OLD_YOSSI_FILTERED']} new_yossi={row['NEW_YOSSI_FILTERED']} "
          f"yossi_violations={row['YOSSI_VIOLATIONS']}")
    return row


def save_row(row: dict):
    # one file per sample: parallel workers each own their own file, so summary.tsv can be
    # assembled afterwards without any of them writing to the same path
    os.makedirs(config.COMPARISON_DIR, exist_ok=True)
    with open(config.comparison_row_path(row["SAMPLE"]), "w") as f:
        json.dump(row, f, indent=2)


def summarize(samples=None) -> bool:
    # merges whatever per-sample rows exist into summary.tsv and reports the verdict
    samples = samples or config.SAMPLES
    rows = []
    for sample in samples:
        fp = config.comparison_row_path(sample)
        if os.path.exists(fp):
            with open(fp) as f:
                rows.append(json.load(f))
    if not rows:
        print("no per-sample comparison rows found -- nothing to summarize")
        return False

    os.makedirs(config.COMPARISON_DIR, exist_ok=True)
    summary_fp = os.path.join(config.COMPARISON_DIR, "summary.tsv")
    pd.DataFrame(rows).to_csv(summary_fp, sep="\t", index=False)
    print(f"\nsummary ({len(rows)}/{len(samples)} samples) -> {summary_fp}")

    missing = [s for s in samples if not os.path.exists(config.comparison_row_path(s))]
    subset_holds = all(r["RAW_SUBSET_OK"] and r["YOSSI_SUBSET_OK"] for r in rows)

    if not subset_holds:
        # a real violation outranks everything else: the subset claim is disproven
        print("OVERALL: FAIL -- see *.raw_violations.tsv / *.yossi_violations.tsv "
              "for the offending loci")
    elif missing:
        # never report PASS on a partial run -- an unfinished sample is not a passing one
        print(f"OVERALL: INCOMPLETE -- the subset claim holds for the {len(rows)} sample(s) that "
              f"finished, but {len(missing)} did not produce results: {', '.join(missing)}")
    else:
        print("OVERALL: PASS -- newest MSMuTect's mutations are a subset of the old version's "
              "for all {} samples".format(len(rows)))
    return subset_holds and not missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--summarize-only", action="store_true",
                         help="skip comparing and just rebuild summary.tsv from existing per-sample rows")
    args = parser.parse_args()

    if not args.summarize_only:
        for sample in args.samples:
            save_row(compare_one(sample))
    summarize(args.samples)


if __name__ == "__main__":
    main()

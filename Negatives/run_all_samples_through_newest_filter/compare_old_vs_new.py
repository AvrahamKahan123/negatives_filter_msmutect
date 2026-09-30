import argparse
import os

import pandas as pd

import config
from results_postprocessing.NoisyLocusDB import NoisyLocusDB
from results_postprocessing.analyze_full_tcga_mut_set import yossi_filter

# Step 4: verify that the newest MSMuTect marks as mutations a SUBSET of what the old
# version did (this is the actual goal stated at the top of Plan.md). Two comparisons are
# made per sample, both keyed on locus = (CHROMOSOME, START, END, PATTERN):
#
#   raw:   new MSMuTect's own calls (--from_file rerun's *.full.mutated_only.mut.tsv,
#          i.e. MSMuTect's internal CALL="M") vs. the old version's calls
#          (*.called.filt.mut.tsv.gz, which is already CALL="M"-only).
#   yossi: the current yossi_filter() applied to each of the above.
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


def new_yossi_filtered(sample: str, noisy_db: NoisyLocusDB) -> pd.DataFrame:
    mutated_only_fp = config.rerun_output_prefix(sample) + ".full.mutated_only.mut.tsv"
    if not os.path.exists(mutated_only_fp):
        raise FileNotFoundError(f"{mutated_only_fp} not found -- run rerun_msmutect_from_file.py first")
    new_mutations_df = pd.read_csv(mutated_only_fp, delimiter="\t")
    filtered = yossi_filter(new_mutations_df.copy(), noisy_db, name=sample)
    os.makedirs(config.NEW_YOSSI_DIR, exist_ok=True)
    filtered.to_csv(config.new_yossi_filtered_path(sample), sep="\t", index=False)
    return new_mutations_df, filtered


def compare_one(sample: str, noisy_db: NoisyLocusDB) -> dict:
    print(f"[{sample}] comparing old vs. new")

    old_mutations_df = pd.read_csv(config.filt_gz_local_path(sample), delimiter="\t")
    old_yossi_fp = config.old_yossi_filtered_path(sample)
    if not os.path.exists(old_yossi_fp):
        raise FileNotFoundError(f"{old_yossi_fp} not found -- run run_yossi_filter_on_old.py first")
    old_yossi_df = pd.read_csv(old_yossi_fp, delimiter="\t")

    new_mutations_df, new_yossi_df = new_yossi_filtered(sample, noisy_db)

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    args = parser.parse_args()

    noisy_db = NoisyLocusDB()
    rows = [compare_one(sample, noisy_db) for sample in args.samples]

    os.makedirs(config.COMPARISON_DIR, exist_ok=True)
    summary_fp = os.path.join(config.COMPARISON_DIR, "summary.tsv")
    pd.DataFrame(rows).to_csv(summary_fp, sep="\t", index=False)
    print(f"\nsummary -> {summary_fp}")

    overall_ok = all(r["RAW_SUBSET_OK"] and r["YOSSI_SUBSET_OK"] for r in rows)
    print("OVERALL: " + ("PASS -- newest MSMuTect's mutations are a subset of the old version's for every sample"
                          if overall_ok else
                          "FAIL -- see *.raw_violations.tsv / *.yossi_violations.tsv for the offending loci"))


if __name__ == "__main__":
    main()

import argparse
import csv
import os
from collections import Counter

import pandas as pd

import config
import msmutect_format
import old_mutation_calls
from results_postprocessing.analyze_full_tcga_mut_set import yossi_filter

# Does (old MSMuTect + yossi_filter) reproduce the CURRENT MSMuTect exactly?
#
# The current version folded yossi's post-processing into the caller itself -- same
# thresholds, same rules (src/IndelCalling/CallMutations.py: normal read support 16 for a
# known germline variant else 9, and 90% of normal reads in the first two alleles). So the
# two call sets ought to coincide. This script measures whether they actually do.
#
#   A = yossi_filter(old run's mutation calls)   <- data/filtered_full_old/*.filtered.full.mut.tsv
#   B = current MSMuTect's mutation calls        <- the --from_file rerun's *.mutated_only.mut.tsv
#
# "Perfectly encapsulates" means A == B. A subset test alone would miss the interesting
# half, so both differences are reported:
#   A \ B  loci the old-version+filter keeps that the current version does NOT call
#   B \ A  loci the current version calls that the old-version+filter does NOT keep
#
# --explain additionally looks up what each version called the disagreeing loci, which is
# what tells you WHY they differ rather than just how much.

KEY = ["CHROMOSOME", "START", "END", "PATTERN"]


def locus_keys(df: pd.DataFrame) -> set:
    # stringify so a column read as int64 in one file and object in another still matches
    return {tuple(str(v) for v in row)
            for row in df[KEY].itertuples(index=False, name=None)}


def calls_for_loci(path: str, wanted: set) -> dict:
    """Stream a *.full.mut.tsv[.gz] and return {locus: CALL} for the loci in `wanted`."""
    found = {}
    if not wanted or not os.path.exists(path):
        return found
    with old_mutation_calls.open_maybe_gzip(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        try:
            idx = [header.index(c) for c in KEY]
            call_index = header.index("CALL")
        except ValueError:
            return found
        for row in reader:
            locus = tuple(row[i] for i in idx)
            if locus in wanted:
                found[locus] = row[call_index]
                if len(found) == len(wanted):
                    break
    return found


def write_loci(loci: set, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(sorted(loci), columns=KEY).to_csv(path, sep="\t", index=False)


def run_one(sample: str, noisy_db=None, explain: bool = False) -> dict:
    source = config.filtered_full_old_path(sample)
    if not os.path.exists(source):
        raise FileNotFoundError(f"{source} not found -- run make_filtered_full_old.py first")
    new_fp = config.new_mutated_only_path(sample)
    if not os.path.exists(new_fp):
        raise FileNotFoundError(f"{new_fp} not found -- run rerun_msmutect_from_file.py first")

    print(f"[{sample}] yossi_filter on {os.path.basename(source)}")
    old_calls = pd.read_csv(source, delimiter="\t", low_memory=False)
    old_filtered = yossi_filter(msmutect_format.normalize_for_yossi_filter(old_calls),
                                noisy_db or config.noisy_locus_db(), name=sample)

    os.makedirs(config.ENCAPSULATION_DIR, exist_ok=True)
    old_filtered.to_csv(config.encapsulation_path(sample, "old_plus_yossi"), sep="\t", index=False)

    new_calls = pd.read_csv(new_fp, delimiter="\t", low_memory=False)
    a, b = locus_keys(old_filtered), locus_keys(new_calls)
    only_old, only_new = a - b, b - a

    write_loci(only_old, config.encapsulation_path(sample, "only_in_old_plus_yossi"))
    write_loci(only_new, config.encapsulation_path(sample, "only_in_new"))

    row = {
        "SAMPLE": sample,
        "OLD_CALLS": len(old_calls),
        "OLD_PLUS_YOSSI": len(a),
        "NEW_CALLS": len(b),
        "SHARED": len(a & b),
        "ONLY_IN_OLD_PLUS_YOSSI": len(only_old),
        "ONLY_IN_NEW": len(only_new),
        "EXACT_MATCH": not only_old and not only_new,
    }

    if explain:
        # what does the OTHER version say about each disagreeing locus?
        print(f"[{sample}]   explaining {len(only_old):,} + {len(only_new):,} disagreements "
              f"(scans both full files, a few minutes)")
        new_says = calls_for_loci(config.new_full_mut_path(sample), only_old)
        old_says = calls_for_loci(config.full_gz_local_path(sample), only_new)
        row["ONLY_IN_OLD_PLUS_YOSSI_BY_NEW_CALL"] = dict(
            Counter(new_says.get(k, "<absent>") for k in only_old).most_common())
        row["ONLY_IN_NEW_BY_OLD_CALL"] = dict(
            Counter(old_says.get(k, "<absent>") for k in only_new).most_common())

    verdict = "EXACT MATCH" if row["EXACT_MATCH"] else "DIFFER"
    print(f"[{sample}] {verdict}: old+yossi={row['OLD_PLUS_YOSSI']:,} new={row['NEW_CALLS']:,} "
          f"shared={row['SHARED']:,} | only_old+yossi={row['ONLY_IN_OLD_PLUS_YOSSI']:,} "
          f"only_new={row['ONLY_IN_NEW']:,}")
    if explain:
        for label, key in (("only in old+yossi, current version calls them",
                            "ONLY_IN_OLD_PLUS_YOSSI_BY_NEW_CALL"),
                           ("only in new, old version called them", "ONLY_IN_NEW_BY_OLD_CALL")):
            if row[key]:
                print(f"[{sample}]   {label}: "
                      + ", ".join(f"{k}={v:,}" for k, v in row[key].items()))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--explain", action="store_true",
                         help="also report what each version called the disagreeing loci "
                              "(scans the full files, adds a few minutes per sample)")
    args = parser.parse_args()

    noisy_db = config.noisy_locus_db()
    rows = [run_one(s, noisy_db, explain=args.explain) for s in args.samples]

    os.makedirs(config.ENCAPSULATION_DIR, exist_ok=True)
    summary_fp = os.path.join(config.ENCAPSULATION_DIR, "encapsulation_summary.tsv")
    pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, dict)} for r in rows]).to_csv(
        summary_fp, sep="\t", index=False)
    print(f"\nsummary -> {summary_fp}")

    exact = [r["SAMPLE"] for r in rows if r["EXACT_MATCH"]]
    if len(exact) == len(rows):
        print(f"VERDICT: yossi_filter on the old version reproduces the current version EXACTLY "
              f"for all {len(rows)} sample(s)")
    else:
        total_old = sum(r["ONLY_IN_OLD_PLUS_YOSSI"] for r in rows)
        total_new = sum(r["ONLY_IN_NEW"] for r in rows)
        print(f"VERDICT: NOT an exact encapsulation -- {len(exact)}/{len(rows)} sample(s) match. "
              f"Across all samples {total_old:,} loci are kept only by old+yossi and "
              f"{total_new:,} only by the current version "
              f"(see *.only_in_old_plus_yossi.tsv / *.only_in_new.tsv)")


if __name__ == "__main__":
    main()

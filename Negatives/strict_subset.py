import sys

# Columns that identify a locus. Both the new (*.full.mutated_only.mut.tsv) and the old
# (*.old_full_mutations.tsv) formats carry these, though the rest of the columns differ.
LOCUS_COLUMNS = ("CHROMOSOME", "START", "END")


def chromosome_rank(chromosome: str):
    """Order chromosomes the way the mut.tsv files are sorted: 1..22, then X, Y, M."""
    name = chromosome[3:] if chromosome.lower().startswith("chr") else chromosome
    if name.isdigit():
        return 0, int(name), ""
    return 1, 0, name.upper()


def sort_key(locus):
    """The order both files are sorted by: chromosome, then END (NOT start)."""
    chromosome, _start, end = locus
    return chromosome_rank(chromosome), end


def read_rows(fp: str):
    """Return (header_line, rows), where each row is ((chromosome, start, end), full_line)."""
    handle = open(fp)
    header = handle.readline().rstrip("\n")
    try:
        columns = [header.split("\t").index(column) for column in LOCUS_COLUMNS]
    except ValueError:
        handle.close()
        raise ValueError(f"{fp} is missing one of {LOCUS_COLUMNS}")

    def rows():
        with handle:
            for line in handle:
                if not line.strip():
                    continue
                line = line.rstrip("\n")
                fields = line.split("\t")
                chromosome, start, end = (fields[column] for column in columns)
                yield (chromosome, int(start), int(end)), line

    return header, rows()


def blocks_by_sort_key(rows):
    """Collapse a sorted row stream into (sort_key, set_of_loci) blocks.

    Only the loci are kept, not the lines -- the superset is far too large to hold in memory,
    and a violation has no superset line to report anyway.

    A block holds every locus sharing a (chromosome, end) -- in practice one, but loci with
    the same end and different starts would otherwise be skipped past by the merge below.
    """
    block_key = None
    block = set()
    for locus, _line in rows:
        key = sort_key(locus)
        if key != block_key:
            if block_key is not None:
                yield block_key, block
            block_key, block = key, set()
        block.add(locus)
    if block_key is not None:
        yield block_key, block


def main(should_be_subset_fp: str, should_be_superset_fp: str):
    ## Algorithm ##
    ## based on the fact that both files are sorted
    ## advance one line at a time in subset
    ## advance line by line in superset until either we hit the locus (chrom, start, end) or we pass it (new chrom or end > superset end). note it is sorted by end and not start
    ## if there is a subset/superset violation, print the offending line from the subset file
    # at the end, print the number of offending loci, and percentage of offending loci of the superset file
    _, superset_rows = read_rows(should_be_superset_fp)
    superset = blocks_by_sort_key(superset_rows)
    superset_key, superset_block = next(superset, (None, set()))

    subset_header, subset_rows = read_rows(should_be_subset_fp)
    print(subset_header)
    subset_total = 0
    offending = 0
    previous_locus = None
    for locus, line in subset_rows:
        subset_total += 1
        key = sort_key(locus)
        if previous_locus is not None and key < sort_key(previous_locus):
            raise ValueError(f"{should_be_subset_fp} is not sorted: {locus} comes after {previous_locus}")
        previous_locus = locus

        # advance the superset until it reaches this locus or passes it
        while superset_key is not None and superset_key < key:
            superset_key, superset_block = next(superset, (None, set()))

        if superset_key == key and locus in superset_block:
            continue
        offending += 1
        print(line)

    print(f"offending loci: {offending}", file=sys.stderr)
    print(f"subset file loci: {subset_total} ({100 * offending / subset_total:.4f}% offending)"
          if subset_total else "subset file loci: 0", file=sys.stderr)

    return offending


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        main(
            "/home/avraham/MaruvkaLab/msmutect_postprocessing/Negatives/run_all_samples_through_newest_filter/results/from_file_rerun/TCGA-A6-5661/TCGA-A6-5661.full.mutated_only.mut.tsv",
            "/home/avraham/MaruvkaLab/msmutect_postprocessing/Negatives/run_all_samples_through_newest_filter/results/old_full_mutations/TCGA-A6-5661.old_full_mutations.tsv"
        )

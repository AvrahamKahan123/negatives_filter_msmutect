#!/usr/bin/env python3
"""Simulate a lower-purity tumor by re-mixing the per-locus read histograms of a real one.

Reads an MSMuTect `*.mut.tsv` -- which carries both the normal and the tumor read histogram
for every locus -- and rewrites the tumor histogram as a mixture of normal (wild-type) reads
and mutant reads at a chosen purity. Sweeping purity gives one sensitivity curve per sample:
how many of this sample's mutations would still be called if the tumor were only X% pure.

Columns are resolved BY NAME from the header, so this reads both MSMuTect layouts:

    old:  ... REFERENCE_REPEATS | NORMAL_MOTIF_REPEATS_* ... | "Noisy Locus" | ...
          unused histogram slots padded with 0, alleles 2-4 spelled NORMAL_ALLELES_n
    new:  ... REFERENCE_REPEATS | KNOWN_GERMLINE_VARIANT | NORMAL_MOTIF_REPEATS_* ...
          unused histogram slots padded with NA, alleles spelled NORMAL_ALLELE_n

The two layouts differ by one added column (KNOWN_GERMLINE_VARIANT) and two removed ones
(the "Noisy Locus" pair). Those happen to cancel out for most -- but not all -- of the fixed
offsets the previous version of this file used, so it read a new-format file silently wrong.

Output is a tumor histogram `.hist.tsv` in the shape `msmutect --from_file` expects (NA
padding, TUMOR_-prefixed histogram columns, locus columns up front), matching what
`running_from_file/split_histograms.py` writes. `--from_file` walks the normal and tumor
histogram files in LOCKSTEP, so the matching normal histogram is written alongside it --
same loci, same order, same number of lines.

    python mix_histograms.py sample.full.mut.tsv out/sample
    python mix_histograms.py sample.full.mut.tsv out/sample --purity 0.3 --replicates 5
"""

import argparse
import gzip
import os
import random
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

# MSMuTect reports at most 6 repeat lengths per histogram, and at most 4 called alleles.
MAX_HISTOGRAM_SLOTS = 6
MAX_ALLELE_SLOTS = 4

# Locus columns, in the order MSMuTect's own Locus.header() writes them.
LOCUS_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN", "REFERENCE_SEQUENCE", "REFERENCE_REPEATS"]
GERMLINE_COLUMN = "KNOWN_GERMLINE_VARIANT"

# What an unused slot looks like on the way out. MSMuTect's reader
# (ResultsReaders.NamedColumns.group) stops at the first "NA", so output must be left-packed.
MISSING = "NA"

CONVERTER_HINT = (
    f"Add it with converting_from_old_format_to_new/convert_old_msmutect_file_to_new.py, "
    f"which drops the old 'Noisy Locus' columns and inserts {GERMLINE_COLUMN} after "
    f"REFERENCE_REPEATS.")


class UnreadableMutFile(RuntimeError):
    """The input does not look like an MSMuTect *.mut.tsv this module can read."""


def open_maybe_gzip(path: str):
    # the lab's files are routinely gzipped
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


class HistogramSet:
    def __init__(self, alleles: List[int] = None, repeat_lengths: Dict[int, int] = None):
        # alleles are the alleles called by the algorithm
        # repeat_lengths is of form {MOTIF: MOTIF_SUPPORT}
        self.alleles = alleles
        if repeat_lengths is None:
            self._repeat_lengths = dict()
        else:
            self._repeat_lengths = repeat_lengths
        self._probability_map: np.ndarray = None # is array of form [REPEAT_LENGTH_1 (REPEAT_LENGTH_1_SUPPORT times), REPEAT_LENGTH_2 (REPEAT_LENGTH_2_SUPPORT times), ...]
        self._sorted_repeat_lengths: Dict[int, int] = None

    def homozygous(self):
        if self.alleles is None:
            raise RuntimeError("No alleles set")
        return len(self.alleles) == 1

    def _invalidate(self):
        self._probability_map = None
        self._sorted_repeat_lengths = None

    def remove_repeat_length(self, length: int):
        if length in self._repeat_lengths:
            del self._repeat_lengths[length]
            self._invalidate()

    def num_reads(self) -> int:
        return sum(self._repeat_lengths.values())

    def num_lengths(self) -> int:
        return len(self._repeat_lengths)

    def add_read(self, repeat_length: int):
        if repeat_length in self._repeat_lengths:
            self._repeat_lengths[repeat_length]+=1
        else:
            self._repeat_lengths[repeat_length]=1
        self._invalidate() # is now invalid

    def create_probability_map(self) -> None:
        pmap = np.zeros((self.num_reads()), dtype=np.int32)
        current_index = 0
        for repeat, repeat_support in self._repeat_lengths.items():
            pmap[current_index:current_index+repeat_support] = repeat
            current_index+=repeat_support
        self._probability_map = pmap

    def randomly_select_read(self) -> int:
        if self.num_reads()==0:
            raise RuntimeError("Cannot randomly select from an empty histogram")
        if self._probability_map is None:
            self.create_probability_map()
        random_int = random.randint(0, self.num_reads()-1)
        return self._probability_map[random_int]

    def copy(self):
        return HistogramSet(alleles=self.alleles, repeat_lengths=self._repeat_lengths.copy())

    def create_sorted_repeat_lengths(self):
        sorted_keys = list(sorted(self._repeat_lengths, key=self._repeat_lengths.get, reverse=True))
        ret = {k: self._repeat_lengths[k] for k in sorted_keys}
        return ret


    @property
    def repeat_lengths(self) -> Dict[int, int]:
        if self._sorted_repeat_lengths is None:
            self._sorted_repeat_lengths = self.create_sorted_repeat_lengths()
        return self._sorted_repeat_lengths


class MutFileSchema:
    """Resolves a *.mut.tsv header to column indices, so rows are read by name not position.

    Mirrors src/Entry/ResultsReaders.py::NamedColumns in MSMuTect itself, including its
    "stop at the first unused slot" rule for the 6-wide histogram groups.
    """

    def __init__(self, header_line: str):
        names = [name.strip() for name in header_line.rstrip("\n").split("\t")]
        self.index = {name: i for i, name in enumerate(names)}
        self.has_germline = GERMLINE_COLUMN in self.index

        # The old format pads unused slots with 0, the new one with NA. Only treat 0 as
        # padding on an old-format file: there a genuine repeat length of 0 is
        # indistinguishable from padding, and padding is overwhelmingly the common case.
        self.blanks = {MISSING, "", "."} if self.has_germline else {MISSING, "", ".", "0"}

        missing_locus = [name for name in LOCUS_COLUMNS if name not in self.index]
        if missing_locus:
            raise UnreadableMutFile(
                f"header is missing {', '.join(missing_locus)}; is this really an MSMuTect "
                f"*.mut.tsv? (columns: {', '.join(names[:8])}...)")

        self.locus_indices = [self.index[name] for name in LOCUS_COLUMNS]
        if self.has_germline:
            self.locus_indices.append(self.index[GERMLINE_COLUMN])

        self.normal_repeats = self._group("NORMAL_MOTIF_REPEATS", MAX_HISTOGRAM_SLOTS)
        self.normal_support = self._group("NORMAL_SUPPORTING_READS", MAX_HISTOGRAM_SLOTS)
        self.tumor_repeats = self._group("TUMOR_MOTIF_REPEATS", MAX_HISTOGRAM_SLOTS)
        self.tumor_support = self._group("TUMOR_SUPPORTING_READS", MAX_HISTOGRAM_SLOTS)
        # the old format spells these NORMAL_ALLELE_1, NORMAL_ALLELES_2, NORMAL_ALLELES_3...
        self.normal_alleles = self._group("NORMAL_ALLELE", MAX_ALLELE_SLOTS,
                                          alternates=("NORMAL_ALLELES",))
        for label, group in (("NORMAL_MOTIF_REPEATS", self.normal_repeats),
                             ("NORMAL_SUPPORTING_READS", self.normal_support),
                             ("TUMOR_MOTIF_REPEATS", self.tumor_repeats),
                             ("TUMOR_SUPPORTING_READS", self.tumor_support),
                             ("NORMAL_ALLELE", self.normal_alleles)):
            if not group:
                raise UnreadableMutFile(f"header has no {label}_* columns")

    def _group(self, base: str, count: int, alternates: Tuple[str, ...] = ()) -> List[int]:
        """Indices of base_1..base_count, stopping at the first one not in the header."""
        indices = []
        for k in range(1, count + 1):
            for candidate in (base,) + alternates:
                name = f"{candidate}_{k}"
                if name in self.index:
                    indices.append(self.index[name])
                    break
            else:
                break
        return indices

    def _values(self, fields: List[str], indices: Iterable[int]) -> List[int]:
        """Read a left-packed group, stopping at the first unused slot."""
        values = []
        for i in indices:
            raw = fields[i].strip()
            if raw in self.blanks:
                break
            values.append(int(float(raw)))
        return values

    def _histogram(self, fields: List[str], repeat_indices, support_indices) -> Dict[int, int]:
        repeats = self._values(fields, repeat_indices)
        support = self._values(fields, support_indices)
        # zip pairs positionally and the shorter list wins, so a ragged row drops its tail
        # rather than pairing a repeat length with some other slot's read count
        return dict(zip(repeats, support))

    def histogram_sets(self, line: str) -> Tuple[HistogramSet, HistogramSet]:
        fields = line.rstrip("\n").split("\t")
        normal = HistogramSet(self._values(fields, self.normal_alleles),
                              self._histogram(fields, self.normal_repeats, self.normal_support))
        tumor = HistogramSet(repeat_lengths=self._histogram(fields, self.tumor_repeats,
                                                            self.tumor_support))
        return normal, tumor

    def locus_fields(self, line: str) -> str:
        fields = line.rstrip("\n").split("\t")
        return "\t".join(fields[i] for i in self.locus_indices)

    def output_header(self, prefix: str) -> str:
        columns = list(LOCUS_COLUMNS)
        if self.has_germline:
            columns.append(GERMLINE_COLUMN)
        columns += [f"{prefix}MOTIF_REPEATS_{k}" for k in range(1, MAX_HISTOGRAM_SLOTS + 1)]
        columns += [f"{prefix}SUPPORTING_READS_{k}" for k in range(1, MAX_HISTOGRAM_SLOTS + 1)]
        return "\t".join(columns)


@dataclass
class MixStats:
    """What the sampler had to compromise on, so a silent degradation shows up in the log."""
    loci: int = 0
    lines: int = 0
    no_mutant_reads: int = 0      # tumor holds only the normal allele: nothing to mix in
    empty_normal: int = 0         # nothing to draw from: tumor histogram copied through
    truncated: int = 0            # mix had >6 repeat lengths; only the 6 best-supported fit
    truncated_reads: int = 0

    def report(self) -> str:
        lines = [f"  {self.lines:,} lines over {self.loci:,} loci"]
        if self.no_mutant_reads:
            lines.append(f"  {self.no_mutant_reads:,} loci had no mutant reads to mix in "
                         f"(tumor histogram holds only the normal allele)")
        if self.empty_normal:
            lines.append(f"  {self.empty_normal:,} loci had an empty normal or tumor histogram "
                         f"(tumor copied through unchanged)")
        if self.truncated:
            lines.append(f"  {self.truncated:,} mixed histograms exceeded "
                         f"{MAX_HISTOGRAM_SLOTS} repeat lengths; dropped "
                         f"{self.truncated_reads:,} reads in the least-supported slots")
        return "\n".join(lines)


def probabilities_to_map(probabilities: List[float]) -> List[float]:
    absolute_probabilities = [probabilities[0]]
    for p in probabilities[1:]:
        absolute_probabilities.append(p+absolute_probabilities[-1])
    return absolute_probabilities


def mix_histograms_probabilistically(histogram_sets: List[HistogramSet], relative_probabilities: List[float], read_support: int) -> HistogramSet:
    if abs(sum(relative_probabilities)-1) > 1e-5:
        raise RuntimeError("Probabilities do not sum to 1")
    probability_map = probabilities_to_map(relative_probabilities)
    mixed_set = HistogramSet()
    for _ in range(read_support):
        p = random.random()
        new_read = None
        for hist_set, prob in zip(histogram_sets, probability_map):
            if p < prob:
                new_read = hist_set.randomly_select_read()
                break
        if new_read is None:
            raise RuntimeError(f"Probability {p} did not fit map {probability_map}. Most likely cause is a rounding error")
        else:
            mixed_set.add_read(new_read)
    return mixed_set


def find_less_dominant_normal_allele(normal: HistogramSet) -> int:
    return list(normal.repeat_lengths.keys())[1] # is normal, so can only have 2 alleles


def mix_histograms(normal: HistogramSet, tumor: HistogramSet, purity: float,
                   stats: MixStats = None) -> Optional[HistogramSet]:
    """Redraw the tumor's reads as a normal/mutant mixture at `purity`.

    Returns None when there is nothing to draw from, which the caller passes through as the
    unmodified tumor histogram -- the line-for-line pairing `msmutect --from_file` relies on
    means a locus can never simply be dropped.
    """
    read_support = tumor.num_reads()
    true_purity = purity/2 # we divide purity by 2 because we assume the mutation is heterozygous
    if normal.num_reads() == 0 or read_support == 0:
        if stats is not None:
            stats.empty_normal += 1
        return None

    # a normal called heterozygous needs a second repeat length to stand in for the allele
    # the mutation replaced; an unfiltered *.full.mut.tsv has loci without one
    heterozygous = not normal.homozygous() and normal.num_lengths() >= 2
    if heterozygous:
        allele_replaced_with_mutation = find_less_dominant_normal_allele(normal)
        normal_without_less_dominant = normal.copy()
        normal_without_less_dominant.remove_repeat_length(allele_replaced_with_mutation)
        tumor_with_normal_allele_removed = tumor.copy()
        tumor_with_normal_allele_removed.remove_repeat_length(allele_replaced_with_mutation)
        sources = [normal, normal_without_less_dominant, tumor_with_normal_allele_removed]
        weights = [1-purity, true_purity, true_purity]
    else:
        normal_allele = normal.alleles[0] if normal.alleles else next(iter(normal.repeat_lengths))
        tumor_with_normal_allele_removed = tumor.copy()
        tumor_with_normal_allele_removed.remove_repeat_length(normal_allele)
        sources = [normal, tumor_with_normal_allele_removed]
        weights = [1-true_purity, true_purity]

    # every source has to be drawable. An exhausted mutant pool means this locus carries no
    # mutation to dilute, so the mixture collapses onto the normal histogram -- which is the
    # right answer for it, and keeps the draw from raising mid-file.
    if any(source.num_reads() == 0 for source in sources):
        if stats is not None:
            stats.no_mutant_reads += 1
        sources, weights = [normal], [1.0]
    return mix_histograms_probabilistically(sources, weights, read_support)


class MemoryEfficientFileWriter:
    """Buffered line writer that only puts the real filename in place on success.

    Output goes to `<path>.partial` and is renamed when the context manager exits cleanly.
    A condor job that is evicted or killed mid-sweep therefore never leaves a truncated
    .hist.tsv behind -- which matters because --skip-existing would otherwise treat a
    half-written file as done, and `--from_file` reads the pair in lockstep, so a short
    tumor file silently mispairs every locus after the cut.
    """

    def __init__(self, output_fp: str, line_threshold: int = 100_000):
        self.output_fp = output_fp
        self.partial_fp = output_fp + ".partial"
        self.output_file = open(self.partial_fp, 'w')
        self.line_threshold = line_threshold
        self.lines_queue = []

    def writeline(self, line: str):
        self.lines_queue.append(line)
        if len(self.lines_queue) > self.line_threshold:
            self.flush_queue()

    def flush_queue(self):
        if self.lines_queue:
            self.output_file.write("\n".join(self.lines_queue)+"\n")
            self.lines_queue = []

    def close(self, succeeded: bool = True):
        if self.output_file.closed:
            return
        if succeeded:
            self.flush_queue()
        self.output_file.close()
        if succeeded:
            os.replace(self.partial_fp, self.output_fp)
        elif os.path.exists(self.partial_fp):
            # never leave a half-written .partial next to the real outputs
            os.remove(self.partial_fp)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *_):
        self.close(succeeded=exc_type is None)

    def __del__(self):
        # reached only if the writer was used without `with`; an interpreter teardown
        # mid-write is exactly the case the .partial name exists to catch
        self.close(succeeded=False)


def fill_in_lists(lst: list, desired_length: int, fill_in_val: object) -> list:
    missing_length = max(desired_length - len(lst), 0) # at least 0
    ret = lst + [fill_in_val for _ in range(missing_length)]
    return ret[:desired_length]


def format_histogram_set(histogram_set: HistogramSet) -> Tuple[str, int]:
    """Render as 6 repeat lengths + 6 read counts, NA-padded. Also returns reads dropped.

    `repeat_lengths` is support-descending, so when a mixture lands on more than
    MAX_HISTOGRAM_SLOTS lengths the slots that fall off are the least-supported ones.
    """
    lengths = list(histogram_set.repeat_lengths.keys())
    support = list(histogram_set.repeat_lengths.values())
    dropped = sum(support[MAX_HISTOGRAM_SLOTS:])
    lengths = fill_in_lists(lengths, MAX_HISTOGRAM_SLOTS, MISSING)
    support = fill_in_lists(support, MAX_HISTOGRAM_SLOTS, MISSING)
    return "\t".join([str(x) for x in lengths+support]), dropped # append the lists


def create_mixed_histograms_file(input_fp: str, output_fp: str, purity: float,
                                 number_of_entries: int = 1,
                                 schema: MutFileSchema = None) -> MixStats:
    """Write the purity-mixed TUMOR histogram for every locus of `input_fp`."""
    stats = MixStats()
    with open_maybe_gzip(input_fp) as opened_mut_file, \
            MemoryEfficientFileWriter(output_fp) as output_file:
        header_line = opened_mut_file.readline()
        if schema is None:
            schema = MutFileSchema(header_line)
        output_file.writeline(schema.output_header("TUMOR_"))
        for line in opened_mut_file:
            if not line.strip():
                continue
            stats.loci += 1
            normal_set, tumor_set = schema.histogram_sets(line)
            line_start = schema.locus_fields(line) + "\t"
            for _ in range(number_of_entries):
                new_tumor_set = mix_histograms(normal_set, tumor_set, purity, stats)
                if new_tumor_set is None:
                    new_tumor_set = tumor_set
                formatted_tumor_set, dropped = format_histogram_set(new_tumor_set)
                if dropped:
                    stats.truncated += 1
                    stats.truncated_reads += dropped
                stats.lines += 1
                output_file.writeline(line_start+formatted_tumor_set)
    return stats


def create_normal_histogram_file(input_fp: str, output_fp: str, number_of_entries: int = 1,
                                 schema: MutFileSchema = None) -> int:
    """Write the companion NORMAL histogram: same loci, same order, same line count.

    `msmutect --from_file` reads the two histogram files in lockstep, so the normal side has
    to be replicated exactly as the tumor side is. It is re-emitted from the parsed histogram
    rather than copied verbatim, so an old-format file's 0-padding comes out as NA.
    """
    lines = 0
    with open_maybe_gzip(input_fp) as opened_mut_file, \
            MemoryEfficientFileWriter(output_fp) as output_file:
        header_line = opened_mut_file.readline()
        if schema is None:
            schema = MutFileSchema(header_line)
        output_file.writeline(schema.output_header("NORMAL_"))
        for line in opened_mut_file:
            if not line.strip():
                continue
            normal_set, _ = schema.histogram_sets(line)
            formatted, _ = format_histogram_set(normal_set)
            line_start = schema.locus_fields(line) + "\t"
            for _ in range(number_of_entries):
                output_file.writeline(line_start+formatted)
                lines += 1
    return lines


def read_schema(input_fp: str) -> MutFileSchema:
    with open_maybe_gzip(input_fp) as handle:
        return MutFileSchema(handle.readline())


def purity_steps(start: float = 0.1, stop: float = 1.0, step: float = 0.05) -> List[float]:
    return [round(x, 2) for x in np.arange(start, stop+step/2, step)]


def create_randomized_purity_files(input_fp: str, output_prefix: str,
                                   number_of_entries_per_purity: int = 1,
                                   purities: List[float] = None,
                                   write_normal: bool = True,
                                   skip_existing: bool = False) -> None:
    """One mixed tumor histogram per purity, plus the one normal histogram they all pair with.

    `skip_existing` leaves finished files alone, so a condor job that died halfway (or a
    sample that failed out of a batch) can be resubmitted without redoing the whole sweep.
    """
    schema = read_schema(input_fp)
    if write_normal and not schema.has_germline:
        print(f"WARNING: {input_fp} has no {GERMLINE_COLUMN} column, so neither will the "
              f"histograms written here, and `msmutect --from_file` will refuse the normal "
              f"one.\n         {CONVERTER_HINT}", file=sys.stderr)

    output_dir = os.path.dirname(os.path.abspath(output_prefix))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    if write_normal:
        normal_fp = f"{output_prefix}_{number_of_entries_per_purity}x.normal.hist.tsv"
        if skip_existing and os.path.exists(normal_fp):
            print(f"normal histogram -> {normal_fp} (exists, skipped)")
        else:
            lines = create_normal_histogram_file(input_fp, normal_fp,
                                                 number_of_entries_per_purity, schema)
            print(f"normal histogram -> {normal_fp} ({lines:,} lines)")

    for purity in (purity_steps() if purities is None else purities):
        output_fp = f"{output_prefix}_{number_of_entries_per_purity}x_{purity}purity.hist.tsv"
        if skip_existing and os.path.exists(output_fp):
            print(f"purity {purity} -> {output_fp} (exists, skipped)")
            continue
        stats = create_mixed_histograms_file(input_fp, output_fp, purity,
                                             number_of_entries_per_purity, schema)
        print(f"purity {purity} -> {output_fp}")
        print(stats.report(), flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="MSMuTect *.mut.tsv (or .gz), old or new format")
    parser.add_argument("output_prefix",
                        help="written as <prefix>_<N>x_<purity>purity.hist.tsv")
    parser.add_argument("--purity", type=float, action="append", dest="purities",
                        help="a single purity to simulate; repeatable. "
                             "Default: 0.1 to 1.0 in steps of 0.05")
    parser.add_argument("--replicates", type=int, default=1,
                        help="re-sampled tumor histograms per locus (default 1). Above 1 the "
                             "normal histogram is replicated to match")
    parser.add_argument("--no-normal-histogram", action="store_true",
                        help="skip the companion normal histogram (you already have one that "
                             "lines up line-for-line)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="leave output files that already exist alone, so an interrupted "
                             "sweep can be resubmitted")
    parser.add_argument("--seed", type=int, help="seed the RNG for a reproducible mixture")
    args = parser.parse_args(argv)

    if args.seed is not None:
        random.seed(args.seed)

    start = time.time()
    try:
        create_randomized_purity_files(args.input, args.output_prefix, args.replicates,
                                       args.purities, not args.no_normal_histogram,
                                       args.skip_existing)
    except (UnreadableMutFile, FileNotFoundError) as e:
        sys.exit(f"ERROR: {e}")
    print(f"done in {time.time()-start:.1f}s")


if __name__ == '__main__':
    main()

# complex_simulation

Simulated **positives**: take a real MSI sample that MSMuTect already called, and ask *how
many of its mutations survive if the tumor had only been X% pure?*

The idea is to rebuild each tumor histogram read-by-read as a mixture of normal reads and
mutant reads, at a controlled purity, write the result out as a histogram file, re-run
MSMuTect on it, and count how many loci still get called. Sweeping purity from 0.1 to 1.0
gives one sensitivity curve per sample — the plots at the end of this directory.

The pipeline is not a single driver script; it is a chain of steps run by hand:

```
  mix_histograms.py          make the purity-mixed tumor histograms   (the real work)
  mix_histograms.sub         + the normal histogram they pair with    5 condor jobs
         |
  msmutect_from_file.sub     call mutations on each simulated pair    95 condor jobs
         |                   msmutect writes .full.call_counts.tsv itself
         |
  results.txt                per-file mutation counts     STILL COLLATED BY HAND
         |
  save_results.py            parse results.txt into a dict
  create_percentages_file.py normalize counts -> percentages table
         |
  plot_lines.py / make_final_graph.py    the curves
```

## `mix_histograms.py` — the core

Reads any MSMuTect `*.mut.tsv` (which carries *both* the normal and the tumor read
histogram for every locus) and emits a `*.hist.tsv` whose tumor histogram has been
re-synthesized at a chosen purity, plus the normal histogram to pair it with.

Reads old- and new-format files alike, and its output goes straight into the current
`msmutect --from_file`.

### What it reads — both formats

`MutFileSchema` resolves the header to column indices **by name**, so the same code reads
the old and the current MSMuTect layouts. It mirrors `NamedColumns` in MSMuTect's own
`src/Entry/ResultsReaders.py`, including its "stop at the first unused slot" rule.

| | old format | current format |
|---|---|---|
| germline column | — | `KNOWN_GERMLINE_VARIANT`, after `REFERENCE_REPEATS` |
| noisy-locus columns | `Noisy Locus` ×2 | dropped |
| unused histogram slot | `0` | `NA` |
| allele columns | `NORMAL_ALLELE_1`, `NORMAL_ALLELES_2..4` | `NORMAL_ALLELE_1..4` |

Per locus it pulls the normal histogram, the normal's **called alleles**, and the tumor
histogram. Padding is format-specific: `0` only counts as padding on an old-format file,
where a genuine repeat length of 0 is indistinguishable from it (there are a couple of real
`0`-length rows per 100k in the TCGA files). Repeat lengths and read counts are paired
positionally and the shorter list wins, so a ragged row drops its tail rather than pairing
a repeat length with another slot's read count.

The old version of this file used **fixed indices**. The added germline column and the two
removed noisy-locus columns cancel out for most of them, which is why a new-format file
appeared to work while the normal histogram was read one column off.

### `HistogramSet`

A locus histogram, `{repeat_length: supporting_reads}`, plus the list of called alleles.

* `create_probability_map()` flattens the histogram into one entry per read
  (`{3: 2, 4: 1}` → `[3, 3, 4]`), so `randomly_select_read()` is a uniform draw over reads
  — i.e. drawing proportionally to read support. The map is invalidated on every
  `add_read`, and lazily rebuilt.
* `repeat_lengths` returns the histogram **sorted by support, descending**. The
  heterozygous logic below, and the 6-slot truncation, both depend on that ordering.
* `homozygous()` is `len(alleles) == 1` — "MSMuTect called one allele here".

### The mixing model

`mix_histograms(normal, tumor, purity)` keeps the original tumor read depth and redraws
that many reads. `true_purity = purity / 2`, because the simulated mutation is assumed
**heterozygous** — only one of the tumor cell's two copies carries it.

*Normal homozygous* — the simple case. Draw each read from:

| source | weight |
|--------|--------|
| the normal histogram (wild-type) | `1 - purity/2` |
| the tumor histogram **with the normal allele removed** | `purity/2` |

Removing the normal allele from the tumor histogram is what makes the second source a
pure "mutant reads only" pool, instead of re-contributing wild-type reads.

*Normal heterozygous* — the less-dominant normal allele (second key of the
support-sorted histogram) is treated as the one the mutation replaced. Three sources:
the full normal (`1 - purity`), the normal without that allele (`purity/2`), and the
tumor without that allele (`purity/2`).

Sampling is `mix_histograms_probabilistically`: cumulative weights via
`probabilities_to_map`, one `random.random()` per read, pick the first bucket it falls in.
It refuses weights that do not sum to 1 within `1e-5`.

**Degenerate loci.** An unfiltered `*.full.mut.tsv` contains loci the model has no answer
for, and a locus can never simply be dropped (see lockstep, below). So:

| situation | what happens |
|-----------|--------------|
| normal called heterozygous but has <2 repeat lengths | treated as homozygous |
| mutant pool empty (tumor holds only the normal allele) | mixture collapses to the normal — this locus carries no mutation to dilute |
| normal or tumor histogram empty | tumor histogram copied through unchanged |

Each is counted and printed in the per-file summary rather than passing silently. Feeding
the old code an unfiltered file raised `Cannot randomly select from an empty histogram`
partway through instead.

### What it writes

A tumor histogram `.hist.tsv` in the shape `msmutect --from_file` expects — the same shape
`running_from_file/split_histograms.py` writes:

```
CHROMOSOME START END PATTERN REFERENCE_SEQUENCE REFERENCE_REPEATS [KNOWN_GERMLINE_VARIANT]
TUMOR_MOTIF_REPEATS_1..6  TUMOR_SUPPORTING_READS_1..6
```

`NA`-padded and left-packed, because MSMuTect's reader stops at the first `NA`. The
germline column is carried through when the input has one. (The old code wrote 18 columns
headed with `NORMAL_`-prefixed names over tumor data, and a truncated header line.)

**The companion normal histogram.** `--from_file` walks the normal and tumor histogram
files **in lockstep**, line N against line N — it never matches on locus. So
`create_randomized_purity_files` also writes the normal histogram: same loci, same order,
same line count, replicated to match `--replicates`. It is re-emitted from the parsed
histogram rather than copied, so an old-format file's `0`-padding comes out as `NA`.

Written once per sweep as `{prefix}_{n}x.normal.hist.tsv`; the tumor files keep their old
names, `{prefix}_{n}x_{purity}purity.hist.tsv`, so `results.txt` parsing still works.

Output is buffered 100k lines at a time by `MemoryEfficientFileWriter`, now a context
manager (it used to rely on `__del__` to flush the tail).

### Running it

```bash
python mix_histograms.py INPUT.mut.tsv OUTPUT_PREFIX           # full 0.1–1.0 sweep
python mix_histograms.py INPUT.mut.tsv OUT --purity 0.3 --replicates 5 --seed 1
```

`.gz` input is read directly.

| flag | effect |
|------|--------|
| `--purity P` | one purity; repeatable. Default: 0.1 → 1.0 step 0.05 |
| `--replicates N` | re-sampled tumor histograms per locus (default 1); the normal is replicated to match |
| `--no-normal-histogram` | skip it — you already have one that lines up line-for-line |
| `--skip-existing` | leave finished files alone, so an interrupted sweep can be resumed |
| `--seed N` | reproducible mixture |

Each file prints a summary: lines, loci, and every compromise the sampler made.

Every file is written to `<name>.partial` and renamed only once it is complete, so an
evicted job never leaves a truncated histogram behind for `--skip-existing` to mistake for
finished work. That matters more than usual here: `--from_file` pairs the two files by line
number, so a short tumor file would mispair every locus after the cut rather than fail.

## Running the sweep on the cluster (`mix_histograms.sub`)

One condor job per sample, each writing all 19 purities plus the normal histogram they
share. The 5 samples are listed inline in the submit file — there is no job list to
regenerate.

```bash
conda activate genomics          # mix_histograms.py needs numpy
./submit_mix.sh
```

| file | role |
|------|------|
| `submit_mix.sh` | **submit through this**, not `condor_submit` directly — see below |
| `mix_histograms.sub` | the submit file; `queue sample in (...)` over the 5 TCGA samples |
| `run_mix_job.sh` | the executable condor runs per sample, mirroring `run_msmutect_over_whole_tcga_again/run_job.sh` |

**Why the wrapper.** Condor does not create the directory for `log`/`output`/`error`, and a
job whose log files cannot be opened goes straight to Hold — all of them, instantly, before
anything runs. The sibling pipeline never hits this because `make_job_lists.py` creates its
log directory as a side effect of building the job list; a sweep that queues from an inline
list has no such step. `submit_mix.sh` creates it, re-asserts `+x` on `run_mix_job.sh`
(the other way a whole batch gets held at once), and submits.

If jobs are held anyway, the hold reason says exactly why:

```bash
condor_q -held -af ClusterId ProcId HoldReason
condor_release <cluster>      # after fixing it; the sweep resumes, finished purities are skipped
```

Input defaults to step 3 of the whole-TCGA rerun —
`<msmutect_results>/<sample>.full.mutated_only.mut.tsv`. Those are current-MSMuTect calls
in the new format, already cut down to `CALL=="M"` rows, which is exactly the set of
positives whose dilution this simulates.

Every path is an environment variable with a default (`getenv = true` passes your shell
through), so a trial run needs no edit:

```bash
OUTPUT_DIR=/storage/.../scratch condor_submit mix_histograms.sub
FORCE=1 condor_submit mix_histograms.sub        # redo purities that already exist
REPLICATES=3 condor_submit mix_histograms.sub
```

`INPUT_DIR`, `INPUT_SUFFIX`, `OUTPUT_DIR`, `LOG_DIR`, `REPLICATES`, `SEED`, `FORCE`,
`PYTHON`. Resubmitting is safe by default — finished purities are skipped.

Budget roughly **half an hour and ~1.2 GB of output per sample** (19 passes over the input,
each writing a file about the input's size). The histograms **cannot be gzipped**:
MSMuTect's `run_from_file` opens them with a plain `open()`.

## Calling mutations on the simulated pairs (`msmutect_from_file.sub`)

One condor job per **(sample, purity)** pair — 5 × 19 = 95, all in parallel. Each runs
`msmutect --from_file` on the mixed tumor histogram against the sample's normal histogram.

```bash
conda activate genomics          # msmutect.sh runs whichever python3 is on PATH
./submit_msmutect.sh
```

| file | role |
|------|------|
| `submit_msmutect.sh` | **submit through this** — builds the job list, creates the log dir, submits |
| `make_msmutect_job_list.sh` | writes `jobs/msmutect_from_file.txt`, one `sample,purity` per line |
| `msmutect_from_file.sub` | `queue sample,purity from jobs/msmutect_from_file.txt` |
| `run_msmutect_job.sh` | the executable condor runs per pair |

Output lands in `<results_dir>/<sample>_1x_<purity>purity.full.mut.tsv`, and MSMuTect
writes `.full.mutated_only.mut.tsv` and **`.full.call_counts.tsv`** next to it — so the
M count per purity is produced for free, no separate counting step.

**The job list is discovered from disk, not regenerated.** The purity in a filename is
whatever python's `f"{round(x, 2)}"` produced — `0.1`, `0.95`, `1.0`. Re-deriving those
strings in bash is a trap: `printf "%g"` turns `1.0` into `1`, and every job for that
purity would look for a file that does not exist. Reading the names back cannot disagree
with the names that exist. It also means a pair only enters the list if **both** its files
are present, so an unfinished mixing run shows up as a short job list rather than as jobs
failing one by one. `make_msmutect_job_list.sh` prints the count per sample and warns when
it is not 19.

Overridable: `MIXED_DIR`, `RESULTS_DIR`, `LOG_DIR`, `REPLICATES`, `MSMUTECT`, `READ_LEVEL`,
`FORCE`, `SAMPLES`. Resubmitting is safe — a pair that already has output is skipped.

```bash
SAMPLES="TCGA-A6-5661" ./submit_msmutect.sh    # one sample's 19 purities
FORCE=1 ./submit_msmutect.sh                   # redo pairs that already have output
```

All 19 jobs of a sample read the **same** normal histogram at once. 95 concurrent readers
is modest, but `/storage` has returned EIO under load before (see `step1_filter.sub`); if
that shows up in the `.err` files, uncomment `max_materialize` in the submit file.

## The downstream scripts

| file | what it does |
|------|--------------|
| `save_results.py` | Parses `results.txt` (`CASE_1x_0.35:  173006`) into `{case: {purity_pct: count}}`. It only `print`s the dict — the output was copy-pasted into `RESULTS_OF_RUN_DEC_3_2025` and into `make_final_graph.py`. Note the `cases` list has `TCGA-AJ-A3BH` twice. Its input, `results.txt`, is still assembled by hand; the sweep now writes a `.full.call_counts.tsv` per pair that nothing yet collates. |
| `create_percentages_file.py` | The other route to the same numbers: reads `../../results/complex_simulation_stats/total_mut_count_table.txt`, divides each count by the sample's hard-coded total mutation count, writes `percentages_table.csv` (one row per purity, one column per case). |
| `plot_lines.py` | Plots `percentages_table.csv` — purity on x, called-mutation fraction on y. The final figure. |
| `make_final_graph.py` | Standalone version of the same plot with the counts and totals pasted inline, so it runs with no input files. x-axis is the raw index 0–18, not purity. |
| `write_condor_script.py` | Emits `tumor_file normal_file output_prefix` lines for the cluster submission, from `../relevant_tsv.txt`. `guide.txt` calls this one "garbage" — it was written against a different directory layout, and its `normal_file` guess (`<case>.called.filt.hist.tsv`) predates `mix_histograms.py` writing the matching normal histogram itself. |
| `read_num_mutations.sh` | A one-line stub (`gunzip -k`), not a working script. |

## Data files

| file | contents |
|------|----------|
| `cases.txt` | The 5 TCGA MSI samples and their total mutation counts. Flagged as coming from the **old** MSMuTect. |
| `results.txt` | Raw per-file called-mutation counts, one line per purity file. Input to `save_results.py`. |
| `RESULTS_OF_RUN_DEC_3_2025` | The parsed dict for the run of 3 Dec 2025, kept as the record of record. |
| `tsv_files.txt` | The full list of generated `*purity.hist.tsv` paths across all samples. |
| `results/` | Output directory, currently empty. |
| `guide.txt` | The one-line-per-file notes this README expands on. |

## Caveats

* **Purity 1.0 is not "the original tumor".** `true_purity = purity/2`, so at purity 1.0
  the simulated tumor is still half wild-type reads — that is what a heterozygous mutation
  in a 100%-pure tumor looks like. The right-hand end of the curve is therefore well below
  1.0 by construction (TCGA-A6-5661: 512,109 of 798,241 = 64% in the Dec 2025 run), not a
  sign that something was lost.
* **Old-format input cannot be run through the current MSMuTect.** `--from_file` refuses a
  normal histogram without `KNOWN_GERMLINE_VARIANT`, and that value cannot be invented from
  the histograms — it comes from the germline DB. `mix_histograms.py` warns and still writes
  the files; convert the input first with
  `converting_from_old_format_to_new/convert_old_msmutect_file_to_new.py`.
* **Mixtures can exceed 6 repeat lengths.** The format has 6 slots, so the least-supported
  ones are dropped and those reads vanish from the simulated depth. It is not rare —
  ~0.8% of loci at purity 0.3 on TCGA-A6-5661, ~1.9% at purity 1.0. The per-file summary
  reports the count; it used to happen silently.
* Two different mutation-count totals are in play — `cases.txt` /
  `make_final_graph.py` (e.g. 798241 for TCGA-A6-5661) versus
  `create_percentages_file.py` (830832). They normalize to slightly different curves.
* Results recorded in `RESULTS_OF_RUN_DEC_3_2025` predate all of this, and came from the
  old MSMuTect on old-format files (`cases.txt` says so too). They are not comparable
  line-for-line with a fresh run.
* The sweep re-reads the input once per purity (19 passes) rather than holding it in
  memory. Fine for the hour-scale genome-wide runs this is used for.

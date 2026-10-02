# run_all_samples_through_newest_filter

Checks whether the **current** version of MSMuTect calls a *subset* of the mutations the
**old** (May 2025 production) version called, on five TCGA samples:

```
TCGA-A6-5661   TCGA-AJ-A3BH   TCGA-AP-A05N   TCGA-FI-A2D4   TCGA-OR-A5LB
```

The method is to take the old run's own per-locus histograms out of its `*.full.mut.tsv.gz`
and feed them back through the current MSMuTect with `--from_file`, so both versions are
compared on *identical input data*. Everything is downloaded to and written under this
directory. See `Plan.md` for the original brief.

Run everything with the `genomics` conda environment:

```bash
conda activate genomics
```

---

## Quick start

```bash
python download_all.py        # fetch all 10 input files (asks for your SSH password ONCE)
python run_parallel.py        # run the whole pipeline, one process per sample (~3-4 h)
```

Results land in `results/comparison/summary.tsv`.

---

## The pipeline

| # | File | What it does |
|---|------|--------------|
| 1 | `download_all.py` | Downloads every input up front |
| 2 | `old_mutation_calls.py` | Extracts the old run's own `CALL=="M"` rows |
| 3 | `run_yossi_filter_on_old.py` | Runs `yossi_filter()` over those old calls |
| 4 | `rerun_msmutect_from_file.py` | Re-runs current MSMuTect via `--from_file` |
| 5 | `compare_old_vs_new.py` | Compares the two call sets, writes the verdict |

`run_all.py` chains 1-5 for one or more samples; `run_parallel.py` runs `run_all.py` once
per sample concurrently.

### `config.py`
Shared settings: the sample list, remote host/paths, every local path, and the
`sys.path` bootstrap that makes `results_postprocessing` importable. **Import this first**
in any new script here — it is what puts the project root on the path.

Also holds `noisy_locus_db()`, a lazily-loaded shared `NoisyLocusDB` (~3.5M loci) so a
process builds it once rather than once per step.

### `ssh_connection.py`
One shared, authenticated SSH connection for all transfers. Without it every `rsync` opens
its own session and prompts for a password again — 10 prompts for 5 samples, spread over
hours. Opens an OpenSSH control master once and points every later `rsync` at its socket.

```bash
python ssh_connection.py --check    # is the shared connection open?
python ssh_connection.py --close    # close it (it also expires on its own after 4h)
```

### `download_samples.py`
Low-level download helpers (`download_filt`, `download_full`). Normally you call
`download_all.py` instead.

### `download_all.py`
Fetches `*.called.filt.mut.tsv.gz` **and** `*.full.mut.tsv.gz` for all five samples, then
closes the SSH connection. Doing this up front matters: a single sample's MSMuTect run takes
hours, so downloading lazily in between would come back to an expired connection — fatal in
a `nohup`'d run with no terminal to re-prompt.

```bash
python download_all.py [--samples S1 S2] [--force] [--keep-ssh]
```

Roughly 0.5 GB per full file, 30-110 MB per filt file.

### `old_mutation_calls.py`
Pulls the old run's own `CALL=="M"` rows out of its full file and caches them.

**This — not `*.called.filt.mut.tsv.gz` — is the baseline for "what the old version
called."** The filt file has already been filtered and is missing real calls (2,442 on
TCGA-A6-5661, 44,829 on TCGA-OR-A5LB). Measuring a subset claim against it would count
those missing calls as new-version violations.

```bash
python old_mutation_calls.py [--samples ...] [--force]
```

### `run_yossi_filter_on_old.py`  (step 2)
Runs the current `yossi_filter()` over the old version's calls, writing
`results/old_yossi_filtered/`.

### `rerun_msmutect_from_file.py`  (step 3)
Carves a pseudo-normal and pseudo-tumor histogram out of each `*.full.mut.tsv.gz` and runs
`msmutect --from_file` on them. Both histograms are produced in a **single streaming pass**
over the gzip, so the multi-GB decompressed file never touches disk.

Adds the `KNOWN_GERMLINE_VARIANT` column, which the old files predate and which the current
`--from_file` run refuses to start without. The value is read **line-for-line** from the
loci file rather than looked up, because a `*.full.mut.tsv` holds exactly one row per locus
in loci-file order (verified: 27,705,407 rows, zero mismatches, zero leftovers). Every row
re-checks that the two files agree on the locus, so a sample that fails to line up aborts
loudly instead of silently mis-annotating the genome.

```bash
python rerun_msmutect_from_file.py [--samples ...] [--force] [--delete-full-gz]
```

### `compare_old_vs_new.py`  (step 4)
Compares both call sets on locus `(CHROMOSOME, START, END, PATTERN)`, in two ways — raw
(each version's unfiltered calls) and yossi (`yossi_filter()` applied to each). Writes a
per-sample JSON row, violation locus lists, and `summary.tsv`.

Verdict is `PASS` only when all samples finished *and* the subset claim held; a partial run
reports `INCOMPLETE`, never `PASS`.

```bash
python compare_old_vs_new.py [--samples ...] [--summarize-only]
```

### `run_all.py`
Steps 1-5 for one or more samples in a single process. Also the unit of work
`run_parallel.py` spawns.

```bash
python run_all.py [--samples ...] [--force] [--skip-summary] [--close-ssh]
```

### `run_parallel.py`
Runs every sample at once, one process each — `--from_file` is single-core by design, so
the way to use the machine is one process per sample. Each worker writes only its own
sample's files, so nothing collides; `summary.tsv` is assembled after they all finish.

Refuses to start if any input is missing, so workers never contend for (or outlive) the SSH
connection.

```bash
python run_parallel.py [--samples ...] [--force] [--max-workers N] [--allow-missing-inputs]
```

Per-sample logs go to `results/logs/<sample>.log`. Expect ~3-4 h wall-clock and ~47 GB peak.

### `msmutect_format.py`
Adapter between the current MSMuTect output format and `yossi_filter()`, which was written
against the older format. Lives here rather than in `results_postprocessing` so the shared
`yossi_filter` keeps behaving identically for every other analysis.

Two incompatibilities, both of which will bite anyone feeding new-format files to old
post-processing code:

1. `FISHER_TEST_P_VALUE` was renamed to **`FISHER_TEST_PVALUE`**. `yossi_filter` asks for
   the old name, gets a `KeyError`, and its bare `except:` turns that into `exit(-1)` with
   the message `FAILED::: <sample>`.
2. Unused allele slots are now written **empty** where the old format wrote `0`.
   `yossi_filter` builds `NORMAL_S` with a chained `+`, so one NaN propagates and the whole
   row drops out silently. On TCGA-A6-5661 that cut 638,341 mutations to **13**, with no
   error — considerably more dangerous than the crash above.

```python
import msmutect_format
df = msmutect_format.normalize_for_yossi_filter(df)   # no-op on old-format frames
```

---

## Standalone tools

### `msmutect_from_file.py`
**The one file to hand to someone who wants to re-call their own `*.full.mut.tsv` with the
current MSMuTect.** Self-contained and stdlib-only — copy it anywhere, it imports nothing
from this project. Does the whole job in one command:

1. carves the pseudo-normal and pseudo-tumor histograms
2. adds `KNOWN_GERMLINE_VARIANT` from the loci file
3. runs `msmutect.sh ... --from_file`
4. prints the resulting call counts

```bash
python msmutect_from_file.py \
    --full-mut-tsv  sample.full.mut.tsv.gz \
    --loci-file     GRCh38...with_germline_variance \
    --output-prefix results/sample_recalled
```

Produces `results/sample_recalled.full.mut.tsv`, `.full.mutated_only.mut.tsv` and
`.full.call_counts.tsv`.

| Flag | Effect |
|------|--------|
| `--msmutect PATH` | where `msmutect.sh` is (default: `$MSMUTECT_PATH`, then `~/MaruvkaLab/MSMuTect/msmutect.sh`, then `PATH`) |
| `--carve-only` | write the histograms and stop, printing the msmutect command to run |
| `--keep-histograms` | keep the carved histograms (several GB genome-wide; deleted on success by default) |
| `--read-level N` | msmutect `-r`: minimum reads to call an allele |
| `--vcf` | also ask msmutect for a VCF |

Verifies `(CHROMOSOME, START, END)` against the loci file on **every** row and exits
non-zero on the first disagreement, rather than mis-annotating the rest of the genome. It
also refuses an `--output-prefix` that would make msmutect overwrite the input, and
resolves `msmutect.sh` *before* the long carve so a bad path fails in a second rather than
an hour.

Needs MSMuTect's own environment (numpy/scipy/pysam) on `PATH`, because `msmutect.sh` runs
whichever `python3` it finds — `conda activate genomics` first.

### `make_filtered_full_old.py`
Writes a mutation-only copy of every full file to `data/filtered_full_old/`, one per input:

```
data/full/TCGA-A6-5661.full.mut.tsv.gz
    -> data/filtered_full_old/TCGA-A6-5661.filtered.full.mut.tsv
```

Keeps the header, so each output is a valid `*.mut.tsv`. ~1 min and ~300 MB per sample;
`--gzip` writes them compressed.

```bash
python make_filtered_full_old.py [--samples ...] [--force] [--gzip]
```

---

## Analyses

### `yossi_filter_on_filtered_full.py`
Asks whether **old MSMuTect + `yossi_filter` == current MSMuTect**. The current version
folded yossi's post-processing into the caller itself (same 9/16 normal-support thresholds,
same 90%-in-first-two-alleles rule), so the two call sets ought to coincide.

Reports both differences, not just a subset test, and `--explain` looks up what each version
called the disagreeing loci.

```bash
python yossi_filter_on_filtered_full.py [--samples ...] [--explain]
```

Writes `results/yossi_encapsulation/`.

### `experiment_strict_rr.py`
Measures what switching `call_verified_locus` to `reversion_to_reference_strict` actually
does, **without** re-running MSMuTect over all 27.7M loci.

Feeds in two groups that were both called `M` by the previous run — the subset *violations*,
and an evenly-spaced *control* sample of mutations the old version also called — so the
repair and its collateral cost are measured on the same run. Histograms are carved for the
selected loci only, streaming the full file and loci file in lockstep so the germline
annotation stays correct on a sparse subset.

```bash
python experiment_strict_rr.py [--sample S] [--control-size N] [--reuse-histograms]
```

Writes `results/strict_rr_experiment/`. Requires `yossi_filter_on_filtered_full.py` to have
produced the violation lists first.

---

## Layout

```
data/
  filt/               downloaded *.called.filt.mut.tsv.gz
  full/               downloaded *.full.mut.tsv.gz
  filtered_full_old/  mutation-only copies of the full files
results/
  old_full_mutations/ the old run's CALL=="M" rows (cached baseline)
  old_yossi_filtered/ step 2 output
  from_file_rerun/    step 3 output, per sample
  new_yossi_filtered/ yossi_filter applied to step 3 output
  comparison/         step 4: summary.tsv, per-sample rows, violation loci
  yossi_encapsulation/ old+yossi vs current comparison
  strict_rr_experiment/ reversion_to_reference_strict measurements
  logs/               per-sample logs from run_parallel.py
```

Nothing under `data/` or `results/` is tracked by git (the repo `.gitignore` keeps only
`*.py`).

---

## Things worth knowing

- **Use the full file's `M` calls as the old baseline, not `called.filt`.** The filt files
  drop a meaningful and *uneven* share of the old version's calls — 0.3% on TCGA-A6-5661 but
  **15%** on TCGA-OR-A5LB.
- **Disk.** A full run peaks around 47 GB with all five samples in parallel. The carved
  histograms (~4 GB/sample) are deleted after each MSMuTect run; `config.ensure_free_disk_gb()`
  aborts below 15 GB free.
- **Runtime.** MSMuTect `--from_file` runs at ~3,100 loci/s on non-mutation loci plus ~6.6 ms
  per mutation locus, so ~3-4 h per sample — dominated by the mutation loci, not the carving
  (~3 min/sample).
- **The subset claim does not hold.** Across the five samples the current version calls
  **4,106** mutations the old version did not. 98.5% of those were `RR` (Reversion to
  Reference) in the old run — `reversion_to_reference()` was rewritten in MSMuTect commit
  `9de7608` (2025-12-03), *after* the May 2025 run. It is a reclassification at the RR
  boundary, not the caller finding new mutations.

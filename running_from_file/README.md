# running_from_file

Re-call an existing MSMuTect `*.full.mut.tsv` with the current MSMuTect, without going back
to the BAMs. The file already holds the per-locus normal and tumor histograms, so they can
be split back out and fed to `msmutect --from_file`.

Three layers, each calling the one below it (see `Plan.md`):

```
  run_on_directory.py      layer 3   a whole directory, N samples in parallel
         |
  run_single_sample.py     layer 2   one sample: a .full.mut.tsv, or a split pair
         |
  split_histograms.py      layer 1   split a .full.mut.tsv into normal + tumor histograms
```

## Before you start: the germline column

Layer 1 **requires** `KNOWN_GERMLINE_VARIANT`, because the current `--from_file` run refuses
to start without it in the normal histogram. Adding that column is a separate job, done by
`converting_from_old_format_to_new/`:

```bash
python -m converting_from_old_format_to_new.convert_old_msmutect_file_to_new \
    --input-dir data/old_format --output-dir data/old_format_wgermline \
    --pattern "*.full.mut.tsv"
```

```python
from converting_from_old_format_to_new.convert_old_msmutect_file_to_new import (
    convert_old_msmutect_file_to_new, convert_directory)
convert_directory("data/old_format", "data/old_format_wgermline")
```

Feed an un-converted file to this package and it stops with a message naming the converter
rather than inventing the value.

Everything needs MSMuTect's own environment (`conda activate genomics`), because
`msmutect.sh` runs whichever `python3` is on `PATH`.

## Layer 1 — `split_histograms.py`

Splits one `*.full.mut.tsv` (or `.gz`) into a pseudo-normal and pseudo-tumor histogram, in a
single pass. Columns are selected **by name**, so it survives columns being added, removed
or reordered.

```bash
python -m running_from_file.split_histograms \
    --full-mut-tsv  sample.full.mut.tsv \
    --output-normal sample.pseudo_normal.hist.tsv \
    --output-tumor  sample.pseudo_tumor.hist.tsv
```

```python
from running_from_file.split_histograms import split_full_mut_tsv, MissingGermlineColumnError
rows = split_full_mut_tsv("sample.full.mut.tsv", "n.hist.tsv", "t.hist.tsv")
```

Runs no MSMuTect. Raises `MissingGermlineColumnError` if the input predates the column.

## Layer 2 — `run_single_sample.py`

Runs `msmutect --from_file` on one sample, from **either** starting point:

```bash
# from a .full.mut.tsv (splits via layer 1 first)
python -m running_from_file.run_single_sample \
    --full-mut-tsv sample.full.mut.tsv --output-prefix results/sample

# from histograms you already have
python -m running_from_file.run_single_sample \
    --normal-histogram n.hist.tsv --tumor-histogram t.hist.tsv --output-prefix results/sample
```

```python
from running_from_file import run_single_sample
run_single_sample.run_from_full_mut_tsv("sample.full.mut.tsv", "results/sample")
run_single_sample.run_from_histograms("n.hist.tsv", "t.hist.tsv", "results/sample")
```

| Flag | Effect |
|------|--------|
| `--msmutect PATH` | where `msmutect.sh` is (default: `$MSMUTECT_PATH`, then `~/MaruvkaLab/MSMuTect/msmutect.sh`, then `PATH`) |
| `--keep-histograms` | keep the split histograms (several GB genome-wide; deleted on success by default) |
| `--read-level N` | msmutect `-r`: minimum reads to call an allele |
| `--vcf` | also ask msmutect for a VCF |

Refuses an `--output-prefix` that would make msmutect overwrite its own input, and resolves
`msmutect.sh` *before* the long split so a bad path fails in a second rather than an hour.
On failure the histograms are left behind so the run can be retried by hand.

## Layer 3 — `run_on_directory.py`

Globs a directory and hands each file to layer 2. Each `--from_file` run is single-core by
design, so the parallelism is across samples.

```bash
python -m running_from_file.run_on_directory \
    --input-dir data/filtered_full_old_wgermline \
    --output-dir results/recalled \
    --workers 6
```

| Flag | Effect |
|------|--------|
| `--pattern` | glob(s) to match (default `*.full.mut.tsv` and `*.full.mut.tsv.gz`) |
| `--workers N` | concurrent msmutect runs (default 6) |
| `--force` | re-run even if the output already exists (otherwise skipped) |

A file that fails — including one that still needs conversion — is reported and skipped;
it does not take down the other workers. Exit status is non-zero if anything failed.

## Callers on top of the pyramid

**`run_samples_through_msmutect.py`** — the GIAB sample set, a thin wrapper over layer 3
with its paths filled in. For same-sample GIAB pairs the calls are false positives, so this
measures the FPR.

**`rerun_msmutect_from_file.py`** — step 3 of the Negatives pipeline
(`Negatives/run_all_samples_through_newest_filter/`). This one is an **adapter, not part of
the pyramid**: its inputs are the raw downloaded `*.full.mut.tsv.gz`, which predate
`KNOWN_GERMLINE_VARIANT`, and converting them first would mean an extra full read+write of a
~5 GB file per sample. So it keeps its own carve, which adds the column in the *same* pass
that splits the histograms — reading it line-for-line from the annotated loci file, with a
per-row check that the two files agree on the locus — and then hands off to layer 2 for the
msmutect invocation. It is imported by `run_all.py` and `experiment_strict_rr.py`.

If you are not on that pipeline, use the pyramid instead.

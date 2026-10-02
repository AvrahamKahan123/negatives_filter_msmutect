# run_msmutect_over_whole_tcga_again

Re-call every TCGA sample with the current MSMuTect, using the histograms already inside the
May-2025 run's output — no BAMs, no realignment. Three condor steps, one job per sample:

```
  <input_dir>/TCGA-W5-AA2K.full.mut.tsv.gz
    1. filter_to_mutations     keep only CALL=="M"        -> <work_root>/filtered/
    2. add_germline            add KNOWN_GERMLINE_VARIANT -> <work_root>/filtered_wgermline/
    3. run_msmutect_from_file  msmutect --from_file       -> <work_root>/msmutect_results/
```

Because step 1 keeps only the loci the **old** run called `M`, step 3 re-calls exactly that
set — which is what makes old-vs-new a subset test.

Steps 2 and 3 do no real work themselves: step 2 calls
`converting_from_old_format_to_new`, step 3 calls the `running_from_file` pyramid. This
directory supplies only the cluster paths and the condor plumbing.

## Setup

Everything configurable is in **`config.json`** — no paths are hardcoded in the scripts or
submit files. Check these before the first submit:

| Key | What |
|-----|------|
| `input_dir` | the May-2025 run's `*.full.mut.tsv.gz` |
| `work_root` | where all three output directories go |
| `msmutect_path` | **`msmutect.sh` on the cluster** — the one value most likely to be wrong |
| `repo_root` | `null` = parent of this directory; correct when jobs run from the checkout |
| `gzip_intermediates` | `true` — keeps `filtered/` and `filtered_wgermline/` ~4x smaller |
| `delete_filtered_after_germline` | `false`; set `true` to drop `filtered/` once converted |

Any key can be overridden at submit time by the upper-cased environment variable, since the
submit files use `getenv = true`:

```bash
WORK_ROOT=/storage/bfe_maruvka/avrahamk/scratch_run condor_submit step1_filter.sub
```

Check what the jobs will see:

```bash
python -m run_msmutect_over_whole_tcga_again.config
```

**Activate MSMuTect's environment before submitting.** `getenv = true` passes your shell to
the jobs, and `msmutect.sh` runs whichever `python3` is on `PATH` — so step 3 fails with
`ModuleNotFoundError: pysam` if you submit from a bare shell.

## Running

Each step is: build the list, submit, wait.

```bash
python -m run_msmutect_over_whole_tcga_again.make_job_lists --step 1
condor_submit step1_filter.sub
# ... wait ...
python -m run_msmutect_over_whole_tcga_again.make_job_lists --step 2
condor_submit step2_germline.sub
# ... wait ...
python -m run_msmutect_over_whole_tcga_again.make_job_lists --step 3
condor_submit step3_msmutect.sub
```

`make_job_lists` writes one sample per line to `jobs/`, which the submit files consume with
`queue sample from jobs/<step>.txt`.

**It lists only the work that remains.** A sample is left out if that step's output already
exists, and also if its *input* doesn't (so a list never contains jobs that must fail). So
after a partial batch, just re-run `make_job_lists` and re-submit — repeat until it reports
0 jobs. Nothing is recomputed. `--include-done` lists everything regardless; `FORCE=1
condor_submit ...` makes the jobs themselves redo finished samples.

Logs go to `<work_root>/logs/step<N>/<sample>.{log,out,err}`.

## Resources

| Step | Memory | Why |
|------|--------|-----|
| 1 filter | 1 GB | streams the input; flat regardless of file size |
| 2 germline | **3 GB** | loads NoisyLocusDB (~3.5M loci in a python set) **per job** |
| 3 msmutect | 2 GB | `--from_file` streams and spills to disk |

Step 2's memory is the one that bites — at 1 GB those jobs are killed. All three are
single-core; `--from_file` is single-core by design, so the parallelism is across samples.

Expect step 3 to take roughly half an hour to a couple of hours per sample, scaling with
that sample's mutation count (~6.6 ms per mutation locus).

## Safety properties

- **Atomic writes.** Steps 1 and 2 write `.partial` and rename only on success, so a job
  evicted mid-write cannot leave a truncated file that the skip-existing check would then
  treat as done.
- **Resumable.** Every step skips a sample whose output exists, so re-submitting a list is
  always safe.
- **Fails loudly.** Missing input, an un-converted file, or a bad `msmutect_path` exit
  non-zero with a message naming the cause, so condor marks the job failed rather than
  recording a silent success.
- **Columns by name.** `CALL` is located from the header, not by index — the index shifts
  between format versions (step 2 drops two columns and adds one).
- **Unreadable inputs are reported, not dropped.** `os.path.exists()` returns False both for
  a file that is absent and for one the filesystem cannot `stat` (an `Input/output error`
  from `/storage`, a dead mount, permissions). Those mean opposite things, so the job lists
  check the errno instead: *missing* is normal sequencing, *unreadable* is a storage fault
  and gets its own loud warning plus `jobs/step<N>_*_unreadable.txt`. Without that
  distinction a batch of EIO files would silently vanish from the final results.

### "N not ready" vs "EXIST BUT CANNOT BE READ"

`not ready` means the previous step hasn't produced that sample's input yet — expected while
a batch is draining; re-run `make_job_lists` later and the number shrinks.

The `*** EXIST BUT CANNOT BE READ ***` block is different and worth acting on: the source
file is listed by the directory but the filesystem errors on it. Those samples are excluded
from the job list and recorded in `jobs/*_unreadable.txt`. Re-run `make_job_lists` once
storage recovers to pick them up.

## Files

| File | |
|------|---|
| `config.json` / `config.py` | all cluster paths; `config.py` also puts `repo_root` on `sys.path` |
| `filter_to_mutations.py` | step 1 |
| `add_germline.py` | step 2 (wraps `converting_from_old_format_to_new`) |
| `run_msmutect_from_file.py` | step 3 (wraps `running_from_file`) |
| `make_job_lists.py` | writes `jobs/*.txt` for `queue ... from` |
| `run_job.sh` | the condor executable: `run_job.sh <module> <sample>` |
| `step{1,2,3}_*.sub` | the submit files |

`run_job.sh` derives `PYTHONPATH` from its own location, so the checkout can live anywhere
on the cluster without editing anything.

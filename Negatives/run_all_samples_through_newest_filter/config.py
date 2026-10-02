import os
import sys

# --- make results_postprocessing importable, the same way analyze_full_tcga_mut_set.py expects ---
# results_postprocessing/cancer_data.py does `from SamplesDB import ...` (bare, not package-relative),
# so results_postprocessing itself must be on sys.path in addition to the project root.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(THIS_DIR, "..", ".."))
RESULTS_POSTPROCESSING_DIR = os.path.join(PROJECT_ROOT, "results_postprocessing")
for _p in (PROJECT_ROOT, RESULTS_POSTPROCESSING_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# --- samples under test ---
SAMPLES = [
    "TCGA-A6-5661",
    "TCGA-AJ-A3BH",
    "TCGA-AP-A05N",
    "TCGA-FI-A2D4",
    "TCGA-OR-A5LB",
]

# --- remote source (see Negatives/run_all_samples_through_newest_filter/Plan.md) ---
REMOTE_HOST = "tech-ui01.hep.technion.ac.il"
REMOTE_USER = "avrahamk"
REMOTE_RUN_DIR = (
    "/storage/bfe_maruvka/gaiafr/Research_project/WGS_SNVs_indels_analysis_project"
    "/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run"
)
REMOTE_FILT_DIR = f"{REMOTE_RUN_DIR}/MSMuTect_called_mut_filt_fixed"
REMOTE_FULL_DIR = REMOTE_RUN_DIR


def remote_filt_path(sample: str) -> str:
    return f"{REMOTE_FILT_DIR}/{sample}.called.filt.mut.tsv.gz"


def remote_full_path(sample: str) -> str:
    return f"{REMOTE_FULL_DIR}/{sample}.full.mut.tsv.gz"


# --- local layout (everything lives under this directory, per Plan.md) ---
DATA_DIR = os.path.join(THIS_DIR, "data")
FILT_DIR = os.path.join(DATA_DIR, "filt")          # downloaded *.called.filt.mut.tsv.gz
FULL_DIR = os.path.join(DATA_DIR, "full")          # downloaded *.full.mut.tsv.gz (one at a time)
# mutation-only copies of FULL_DIR, made by make_filtered_full_old.py
FILTERED_FULL_OLD_DIR = os.path.join(DATA_DIR, "filtered_full_old")

RESULTS_DIR = os.path.join(THIS_DIR, "results")
# the old run's own CALL=="M" rows, pulled out of its full file (see old_mutation_calls.py).
# This -- not the pre-filtered *.called.filt.mut.tsv.gz -- is the baseline "what the old
# version called", for both yossi_filter and the subset comparison.
OLD_FULL_MUTATIONS_DIR = os.path.join(RESULTS_DIR, "old_full_mutations")
OLD_YOSSI_DIR = os.path.join(RESULTS_DIR, "old_yossi_filtered")     # step 2 output
RERUN_DIR = os.path.join(RESULTS_DIR, "from_file_rerun")            # step 3 output (msmutect --from_file)
NEW_YOSSI_DIR = os.path.join(RESULTS_DIR, "new_yossi_filtered")     # yossi_filter applied to step 3 output
COMPARISON_DIR = os.path.join(RESULTS_DIR, "comparison")            # step 4 output
# does (old MSMuTect + yossi_filter) reproduce the current MSMuTect exactly?
# see yossi_filter_on_filtered_full.py
ENCAPSULATION_DIR = os.path.join(RESULTS_DIR, "yossi_encapsulation")
LOGS_DIR = os.path.join(RESULTS_DIR, "logs")                        # per-sample logs from run_parallel.py

MSMUTECT_PATH = os.path.expanduser("~/MaruvkaLab/MSMuTect/msmutect.sh")

# The phobos loci file annotated with HOMOGENEOUS/HETEROGENEOUS in its last column. A
# *.full.mut.tsv holds one row per locus in exactly this file's order, so the germline
# annotation can be read off line-for-line instead of looked up (see
# rerun_msmutect_from_file.py::carve_histograms).
GERMLINE_LOCI_FILE = ("/home/avraham/MaruvkaLab/msmutect_development/data/"
                      "GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed.with_germline_variance")
# 0-indexed fields in that file: chromosome (with "chr" prefix), start, end, pattern, and
# the HOMOGENEOUS/HETEROGENEOUS annotation (last column)
LOCI_FILE_CHROMOSOME, LOCI_FILE_START, LOCI_FILE_END, LOCI_FILE_PATTERN = 0, 3, 4, 12
LOCI_FILE_GERMLINE = 16

# locus columns carved out for --from_file. KNOWN_GERMLINE_VARIANT is REQUIRED by the
# normal histogram; the old files predate it, so rerun_msmutect_from_file.py reads it
# line-for-line off GERMLINE_LOCI_FILE. Other absent columns are simply skipped.
LOCUS_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN",
                 "REFERENCE_SEQUENCE", "REFERENCE_REPEATS", "KNOWN_GERMLINE_VARIANT"]

# locus key used everywhere to compare/join mutations between old and new runs
LOCUS_KEY_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN"]

# refuse to start a new download/decompression step if free disk space drops below this
MIN_FREE_DISK_GB = 15


def filt_gz_local_path(sample: str) -> str:
    return os.path.join(FILT_DIR, f"{sample}.called.filt.mut.tsv.gz")


def full_gz_local_path(sample: str) -> str:
    return os.path.join(FULL_DIR, f"{sample}.full.mut.tsv.gz")


def rerun_sample_dir(sample: str) -> str:
    return os.path.join(RERUN_DIR, sample)


def rerun_output_prefix(sample: str) -> str:
    return os.path.join(rerun_sample_dir(sample), sample)


def old_full_mutations_path(sample: str) -> str:
    return os.path.join(OLD_FULL_MUTATIONS_DIR, f"{sample}.old_full_mutations.tsv")


def filtered_full_old_path(sample: str) -> str:
    # TCGA-A6-5661.full.mut.tsv[.gz] -> TCGA-A6-5661.filtered.full.mut.tsv
    return os.path.join(FILTERED_FULL_OLD_DIR, f"{sample}.filtered.full.mut.tsv")


def encapsulation_path(sample: str, what: str) -> str:
    return os.path.join(ENCAPSULATION_DIR, f"{sample}.{what}.tsv")


def new_mutated_only_path(sample: str) -> str:
    return rerun_output_prefix(sample) + ".full.mutated_only.mut.tsv"


def new_full_mut_path(sample: str) -> str:
    return rerun_output_prefix(sample) + ".full.mut.tsv"


def old_yossi_filtered_path(sample: str) -> str:
    return os.path.join(OLD_YOSSI_DIR, f"{sample}.old_yossi_filtered.tsv")


def new_yossi_filtered_path(sample: str) -> str:
    return os.path.join(NEW_YOSSI_DIR, f"{sample}.new_yossi_filtered.tsv")


def comparison_row_path(sample: str) -> str:
    # one file per sample so parallel workers never write the same file; summary.tsv is
    # assembled from these afterwards (compare_old_vs_new.summarize)
    return os.path.join(COMPARISON_DIR, f"{sample}.comparison_row.json")


def sample_log_path(sample: str) -> str:
    return os.path.join(LOGS_DIR, f"{sample}.log")


_noisy_db = None


def noisy_locus_db():
    # Shared lazily-loaded NoisyLocusDB (~3.5M loci). yossi_filter needs it for its
    # per-locus normal-support threshold, and both step 2 and step 4 call yossi_filter, so
    # without this each would build its own copy.
    global _noisy_db
    if _noisy_db is None:
        from converting_from_old_format_to_new.NoisyLocusDB import NoisyLocusDB
        print("loading NoisyLocusDB...")
        _noisy_db = NoisyLocusDB()
        print(f"loaded {len(_noisy_db.db)} loci")
    return _noisy_db


def ensure_free_disk_gb(min_gb: float = MIN_FREE_DISK_GB, path: str = THIS_DIR):
    import shutil
    free_gb = shutil.disk_usage(path).free / (1024 ** 3)
    if free_gb < min_gb:
        raise RuntimeError(
            f"Only {free_gb:.1f} GB free at {path} (minimum required: {min_gb} GB). "
            f"Aborting to avoid filling up the disk -- clean up data/ or results/ and retry."
        )

import pandas as pd

import config  # noqa: F401  -- imported for the sys.path setup results_postprocessing needs
from results_postprocessing.enums import COLUMN

# Adapter between the CURRENT MSMuTect output format and yossi_filter(), which was written
# against the older format. Kept here rather than in results_postprocessing so the shared
# yossi_filter keeps behaving exactly as it does for every other analysis.
#
# Two incompatibilities, both found by comparing a from_file rerun against the old files:
#
# 1. FISHER_TEST_P_VALUE was renamed to FISHER_TEST_PVALUE (no underscore). yossi_filter
#    asks for the old name, so it raises KeyError -- which its bare `except:` turns into
#    `exit(-1)` and the one-line message "FAILED::: <sample>".
#
# 2. Unused allele slots are now written EMPTY where the old format wrote 0. yossi_filter
#    computes NORMAL_S with a chained `+`, so a single NaN propagates and NORMAL_S becomes
#    NaN for nearly every row -- both the normal-support and third-and-later tests then
#    silently evaluate False. On TCGA-A6-5661 that cut 638,341 mutations down to 13 with no
#    error at all, which is far more dangerous than the crash above.
#
# Both are format drift, not real signal, so the fix is to restore the old spellings and
# the old "absent allele == 0" convention before handing the frame to yossi_filter.

RENAMED_COLUMNS = {"FISHER_TEST_PVALUE": COLUMN.FISHER_TEST_P_VALUE}

# every per-allele slot: absent used to mean 0, and yossi_filter still reads it that way
ALLELE_SLOT_PREFIXES = ["NORMAL_", "TUMOR_"]
ALLELE_SLOT_BASES = ["MOTIF_REPEATS_", "SUPPORTING_READS_", "ALLELE_", "ALLELES_", "FRACTION_"]
ALLELE_SLOT_RANGE = range(1, 7)


def allele_slot_columns(df: pd.DataFrame):
    wanted = {f"{prefix}{base}{i}"
              for prefix in ALLELE_SLOT_PREFIXES
              for base in ALLELE_SLOT_BASES
              for i in ALLELE_SLOT_RANGE}
    return [c for c in df.columns if c in wanted]


def normalize_for_yossi_filter(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of `df` that yossi_filter can read, whichever format it came from.

    Safe to call on old-format frames: the rename only fires if the new spelling is
    present, and old-format files have no NaNs in the allele slots to fill.
    """
    out = df.copy()

    renames = {old: new for old, new in RENAMED_COLUMNS.items()
               if old in out.columns and new not in out.columns}
    if renames:
        out = out.rename(columns=renames)

    slots = allele_slot_columns(out)
    if slots:
        out[slots] = out[slots].fillna(0)
    return out


def describe_normalization(df: pd.DataFrame) -> str:
    """One-line summary of what normalize_for_yossi_filter would change (for logging)."""
    renames = [old for old, new in RENAMED_COLUMNS.items()
               if old in df.columns and new not in df.columns]
    slots = allele_slot_columns(df)
    filled = int(df[slots].isna().sum().sum()) if slots else 0
    if not renames and not filled:
        return "already in the format yossi_filter expects"
    parts = []
    if renames:
        parts.append("renamed " + ", ".join(f"{r} -> {RENAMED_COLUMNS[r]}" for r in renames))
    if filled:
        parts.append(f"filled {filled:,} empty allele-slot cells with 0")
    return "; ".join(parts)

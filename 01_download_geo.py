"""01_download_geo.py
================================================================================
Download GSE18123, GSE25507, GSE42133 series matrices from NCBI GEO and parse
sample metadata into per-platform phenotype CSVs.

GSE18123 contains samples on TWO Affymetrix platforms (GPL570 and GPL6244).
The pipeline treats these as separate cohorts (training and primary external).
A single sample is uniquely identified by its GSM accession; samples are
assigned to a sub-cohort by the platform field in the series matrix.

Outputs
-------
data/raw/<GSE>_<GPL>_series_matrix.txt.gz   (raw GEO series matrices)
data/raw/<GPL>_platform.txt                 (platform-level probe→gene table)
data/processed/phenotype_<KEY>.csv          (sample_id, asd, sex, age)
data/raw/expr_log2_<KEY>.parquet            (probe-level log2 expression)
"""
from __future__ import annotations
import gzip
import io
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

import config


GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _series_matrix_url(gse: str, platform: str | None = None) -> str:
    # GEO stores series matrices under series/<GSEnnn>/<GSE>/matrix/
    head = gse[:-3] + "nnn"
    fname = f"{gse}_series_matrix.txt.gz"
    if platform is not None:
        fname = f"{gse}-{platform}_series_matrix.txt.gz"
    return f"{GEO_BASE}/{head}/{gse}/matrix/{fname}"


def _platform_table_url(gpl: str) -> str:
    head = gpl[:-3] + "nnn" if len(gpl) > 3 else "GPLnnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{head}/{gpl}/annot/{gpl}.annot.gz"


def download(url: str, dest: Path, retries: int = 3) -> None:
    """Stream a URL to disk; retries on transient errors."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached: {dest.name}")
        return
    for attempt in range(1, retries + 1):
        try:
            print(f"  downloading {url}")
            r = requests.get(url, stream=True, timeout=60)
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 16):
                    f.write(chunk)
            return
        except Exception as e:
            print(f"  attempt {attempt}/{retries} failed: {e}")
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"download failed after {retries} attempts: {url}")


def _open_text(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse a GEO series matrix into (expression, phenotype) DataFrames.

    Returns
    -------
    expr : DataFrame, probes (rows) × samples (cols)
    pheno : DataFrame, one row per sample, columns include
            sample_id, characteristics_*, platform
    """
    meta_rows = {}
    expr_lines = []
    in_table = False
    samples = []
    with _open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                expr_lines.append(line)
                continue
            if line.startswith("!"):
                # Sample-level metadata is repeated as columns
                parts = line.split("\t")
                key = parts[0].lstrip("!")
                values = [v.strip('"') for v in parts[1:]]
                meta_rows.setdefault(key, []).append(values)

    if not expr_lines:
        raise RuntimeError(f"no expression table found in {path}")

    # First line of the expression block is the header: ID_REF \t GSMxxxx ...
    header = expr_lines[0].split("\t")
    samples = [s.strip('"') for s in header[1:]]
    # Body
    body = "\n".join(expr_lines[1:])
    expr = pd.read_csv(io.StringIO(body), sep="\t", header=None,
                        names=header, index_col=0, dtype={0: str},
                        na_values=["", "NA", "null"])
    expr.index = expr.index.map(lambda x: str(x).strip('"'))
    expr.columns = samples
    # The data is samples × probes is more useful for our downstream;
    # but we keep probes × samples here, transpose later in 02_preprocess.
    expr = expr.apply(pd.to_numeric, errors="coerce")

    # Build phenotype frame
    pheno = pd.DataFrame({"sample_id": samples})
    for key, value_lists in meta_rows.items():
        if not value_lists:
            continue
        if len(value_lists[0]) != len(samples):
            continue
        if len(value_lists) == 1:
            pheno[key] = value_lists[0]
        else:
            # multi-line characteristic, e.g. Sample_characteristics_ch1
            for i, vl in enumerate(value_lists):
                pheno[f"{key}_{i + 1}"] = vl

    return expr, pheno


# ---------------------------------------------------------------------------
# Phenotype parsing
# ---------------------------------------------------------------------------

ASD_REGEX = re.compile(r"\b(autism|asd|autistic)\b", flags=re.I)
CTRL_REGEX = re.compile(r"\b(typically developing|control|td|ctl|healthy)\b", flags=re.I)
MALE_REGEX = re.compile(r"\b(male|^m$)\b", flags=re.I)
FEMALE_REGEX = re.compile(r"\b(female|^f$)\b", flags=re.I)
AGE_REGEX = re.compile(r"(\d+(?:\.\d+)?)\s*(?:y|year|yr|mo|month)?", flags=re.I)


def _row_text(row: pd.Series, prefixes=("Sample_characteristics_ch1",
                                          "Sample_source_name_ch1",
                                          "Sample_title",
                                          "Sample_description",)) -> str:
    """Concatenate all candidate annotation fields for a sample row."""
    parts = []
    for col in row.index:
        if any(col.startswith(p) for p in prefixes):
            v = row[col]
            if isinstance(v, str):
                parts.append(v)
    return " | ".join(parts).lower()


def derive_asd_label(row: pd.Series) -> int | None:
    txt = _row_text(row)
    if not txt:
        return None
    # ASD label first so 'autistic' wins over 'autism control' (no such phrasing exists in these series)
    if ASD_REGEX.search(txt):
        return 1
    if CTRL_REGEX.search(txt):
        return 0
    return None


def derive_sex(row: pd.Series) -> str | None:
    txt = _row_text(row)
    if FEMALE_REGEX.search(txt):
        return "F"
    if MALE_REGEX.search(txt):
        return "M"
    return None


def derive_age(row: pd.Series) -> float | None:
    txt = _row_text(row)
    m = re.search(r"age[:\s]+([\d.]+)", txt)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def derive_platform(row: pd.Series) -> str | None:
    for col in row.index:
        if col.lower() == "sample_platform_id":
            v = row[col]
            if isinstance(v, str) and v.startswith("GPL"):
                return v
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_cohort(key: str) -> None:
    info = config.COHORTS[key]
    gse, gpl = info["gse"], info["platform"]
    print(f"[{key}]  gse={gse} platform={gpl} role={info['role']}")

    # GEO sometimes splits a GSE series across platforms with -GPLxxx in the filename
    candidates = [
        config.RAW / f"{gse}-{gpl}_series_matrix.txt.gz",
        config.RAW / f"{gse}_series_matrix.txt.gz",
    ]
    series_path = None
    for c in candidates:
        if c.exists():
            series_path = c
            break
    if series_path is None:
        # Try platform-tagged first (multi-platform series)
        try:
            url = _series_matrix_url(gse, gpl)
            dest = candidates[0]
            download(url, dest)
            series_path = dest
        except Exception:
            url = _series_matrix_url(gse)
            dest = candidates[1]
            download(url, dest)
            series_path = dest

    expr_probes_x_samples, pheno = parse_series_matrix(series_path)

    # If the series matrix mixes platforms, keep only samples on the requested platform
    if "Sample_platform_id" in pheno.columns:
        before = len(pheno)
        keep = pheno["Sample_platform_id"] == gpl
        pheno = pheno[keep].reset_index(drop=True)
        expr_probes_x_samples = expr_probes_x_samples[pheno["sample_id"].tolist()]
        print(f"  platform filter: kept {len(pheno)}/{before} samples on {gpl}")

    # Derive analytic phenotype fields
    pheno["asd"] = pheno.apply(derive_asd_label, axis=1)
    pheno["sex"] = pheno.apply(derive_sex, axis=1)
    pheno["age"] = pheno.apply(derive_age, axis=1)

    n_total = len(pheno)
    n_asd = int((pheno["asd"] == 1).sum())
    n_ctrl = int((pheno["asd"] == 0).sum())
    n_lost = int(pheno["asd"].isna().sum())
    print(f"  derived: n={n_total}, ASD={n_asd}, control={n_ctrl}, unlabelled={n_lost}")

    # Drop samples without a derivable ASD label
    pheno = pheno[pheno["asd"].notna()].copy()
    pheno["asd"] = pheno["asd"].astype(int)
    expr_probes_x_samples = expr_probes_x_samples[pheno["sample_id"].tolist()]

    # Apply analytic restrictions
    if key == "GPL570":
        # Per Methods: training cohort is restricted to males (where sex is annotated)
        mask = pheno["sex"].fillna("M").eq("M")
        n_before = len(pheno)
        pheno = pheno[mask].copy()
        expr_probes_x_samples = expr_probes_x_samples[pheno["sample_id"].tolist()]
        print(f"  male-only restriction: kept {len(pheno)}/{n_before}")

    # Save phenotype
    pheno_min = pheno[["sample_id", "asd", "sex", "age"]].copy()
    pheno_out = config.PROCESSED / f"phenotype_{key}.csv"
    pheno_min.to_csv(pheno_out, index=False)
    print(f"  wrote {pheno_out}")

    # Save log2 expression (probes × samples → samples × probes) as parquet
    expr_samples_x_probes = expr_probes_x_samples.T
    expr_samples_x_probes.index.name = "sample_id"
    expr_out = config.RAW / f"expr_log2_{key}.parquet"
    # Sanity check: GEO series matrices are conventionally on log2 scale already
    # (RMA/fRMA for Affymetrix; neqc for Illumina). Verify by max value heuristic.
    max_val = float(expr_samples_x_probes.values[~pd.isna(expr_samples_x_probes.values)].max())
    if max_val > 30:
        print(f"  WARNING: max expression {max_val:.2f} > 30, suggests raw (non-log) scale; applying log2(x+1)")
        expr_samples_x_probes = (expr_samples_x_probes + 1).apply(lambda c: c.clip(lower=0))
        import numpy as np
        expr_samples_x_probes = expr_samples_x_probes.apply(lambda c: pd.Series(np.log2(c.values + 1), index=c.index))
    else:
        print(f"  log2 scale confirmed (max = {max_val:.2f})")
    expr_samples_x_probes.to_parquet(expr_out)
    print(f"  wrote {expr_out}  shape={expr_samples_x_probes.shape}")


def download_platform_table(gpl: str) -> None:
    dest = config.RAW / f"{gpl}_annot.gz"
    if dest.exists() and dest.stat().st_size > 0:
        return
    try:
        download(_platform_table_url(gpl), dest)
    except Exception as e:
        print(f"  WARNING: could not download {gpl} annotation file: {e}")
        print(f"  Will rely on annotation embedded in the platform soft file (handled in 02_preprocess.py)")


def main():
    for key in config.ALL_KEYS:
        process_cohort(key)

    # Platform-level annotations needed for probe-to-gene mapping in 02_preprocess
    for gpl in sorted({info["platform"] for info in config.COHORTS.values()}):
        download_platform_table(gpl)

    print("\nstage 01 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)

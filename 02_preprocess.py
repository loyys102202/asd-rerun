"""02_preprocess.py
================================================================================
Map probes → genes and aggregate each cohort's expression matrix.

For each gene with multiple mapped probes on a given platform, the probe with
the largest interquartile range across samples in that cohort is retained
(IQR-max-per-gene). This retains the highest-informative-variance probe under
the assumption that low-variance probes are dominated by background noise.

Inputs
------
data/raw/expr_log2_<KEY>.parquet         (sample × probe; from 01)
data/raw/<GPL>_annot.gz                  (platform annotation file from GEO)

Outputs
-------
data/processed/probe_map_<GPL>.csv       (probe, gene)
data/processed/expr_gene_<KEY>.parquet   (sample × gene)
"""
from __future__ import annotations
import gzip
import io
import re
import sys
from pathlib import Path

import pandas as pd

import config
from pipeline import iqr_max_aggregate


# ---------------------------------------------------------------------------
# Platform annotation parsing
# ---------------------------------------------------------------------------
# GEO .annot.gz files are a SOAP-like tab-separated table; the relevant columns
# are 'ID' (probe id) and 'Gene symbol' (or similar). Different platforms use
# slightly different column names, so we probe for the likely candidates.

GENE_COL_CANDIDATES = [
    "Gene symbol",
    "Gene Symbol",
    "gene_assignment",       # GPL6244 native
    "Symbol",                # GPL10558 (Illumina)
    "ILMN_Gene",             # GPL10558 backup
    "GENE_SYMBOL",
]


def _open_annot(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"),
                                 encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_geo_annot(path: Path) -> pd.DataFrame:
    """Parse a GEO platform annotation file into a (probe, gene) table.

    Robust to the various platform-specific header conventions.
    """
    header_line = None
    rows = []
    with _open_annot(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("!") or line.startswith("^"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if header_line is None:
                header_line = line.split("\t")
                continue
            cols = line.split("\t")
            # Pad short rows
            if len(cols) < len(header_line):
                cols += [""] * (len(header_line) - len(cols))
            rows.append(cols[:len(header_line)])

    if header_line is None:
        raise RuntimeError(f"failed to parse {path}")

    df = pd.DataFrame(rows, columns=header_line)
    df.columns = [c.strip() for c in df.columns]

    # Probe id column is conventionally "ID"
    id_col = "ID"
    if id_col not in df.columns:
        for c in df.columns:
            if c.lower() in ("id", "probe_id", "probeid"):
                id_col = c
                break

    # Find gene column
    gene_col = None
    for cand in GENE_COL_CANDIDATES:
        if cand in df.columns:
            gene_col = cand
            break
    if gene_col is None:
        raise RuntimeError(f"could not find gene-symbol column in {path}; "
                            f"saw columns {df.columns.tolist()}")

    out = df[[id_col, gene_col]].rename(columns={id_col: "probe", gene_col: "gene"})
    out["probe"] = out["probe"].astype(str).str.strip()
    out["gene"] = out["gene"].astype(str).str.strip()

    # gene_assignment on GPL6244 has the form 'NM_xxx // SYMBOL // description // ...'
    # The gene symbol is the second '//'-delimited field.
    if gene_col == "gene_assignment":
        out["gene"] = out["gene"].apply(
            lambda s: s.split("//")[1].strip() if "//" in s and len(s.split("//")) >= 2
                      else ""
        )

    # Multi-symbol cells (Illumina sometimes has 'SYMBOL1 /// SYMBOL2'): take first
    out["gene"] = out["gene"].apply(
        lambda s: s.split("///")[0].strip() if isinstance(s, str) else ""
    )

    # Drop probes with no gene symbol
    out = out[out["gene"].astype(bool) & out["gene"].ne("---")]
    out = out.drop_duplicates(subset=["probe"]).reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def process_cohort(key: str, probe_map: pd.DataFrame) -> None:
    expr_path = config.RAW / f"expr_log2_{key}.parquet"
    if not expr_path.exists():
        raise FileNotFoundError(f"missing {expr_path} — run 01_download_geo.py first")
    print(f"[{key}]  loading {expr_path.name}")
    expr = pd.read_parquet(expr_path)
    print(f"  shape (sample × probe) = {expr.shape}")

    # Note: pipeline.iqr_max_aggregate expects samples × probes? No — read carefully:
    # the helper as written takes a DataFrame with samples in ROWS and probes in COLS.
    # Our expr is already samples × probes, so we pass it directly.
    probe_to_gene = probe_map.set_index("probe")["gene"]

    print(f"  aggregating to gene level by IQR-max...")
    expr_gene = iqr_max_aggregate(expr, probe_to_gene)
    print(f"  shape (sample × gene) = {expr_gene.shape}")

    out = config.PROCESSED / f"expr_gene_{key}.parquet"
    expr_gene.to_parquet(out)
    print(f"  wrote {out}")


def main():
    # Resolve unique platforms and load annotations once
    platforms = sorted({info["platform"] for info in config.COHORTS.values()})
    probe_maps: dict[str, pd.DataFrame] = {}
    for gpl in platforms:
        annot = config.RAW / f"{gpl}_annot.gz"
        if not annot.exists():
            raise FileNotFoundError(
                f"missing {annot}; run 01_download_geo.py or download "
                f"{gpl}.annot.gz from "
                f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/.../{gpl}/annot/")
        print(f"[platform]  parsing annotation for {gpl}")
        pm = parse_geo_annot(annot)
        print(f"  {len(pm)} probes mapped to {pm['gene'].nunique()} unique genes")
        pm.to_csv(config.PROCESSED / f"probe_map_{gpl}.csv", index=False)
        probe_maps[gpl] = pm

    for key in config.ALL_KEYS:
        gpl = config.COHORTS[key]["platform"]
        process_cohort(key, probe_maps[gpl])

    print("\nstage 02 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)

"""03_common_genes.py
================================================================================
Find the intersection of genes present on BOTH GPL570 and GPL6244 after
IQR-max probe-to-gene aggregation. This intersection is the LASSO candidate
pool that ensures every panel member is mappable on the primary external
platform.

Inputs
------
data/processed/expr_gene_GPL570.parquet
data/processed/expr_gene_GPL6244.parquet

Outputs
-------
data/processed/common_genes_570_6244.txt
"""
from __future__ import annotations
import sys

import pandas as pd

import config


def main():
    g570 = pd.read_parquet(config.PROCESSED / "expr_gene_GPL570.parquet").columns
    g6244 = pd.read_parquet(config.PROCESSED / "expr_gene_GPL6244.parquet").columns

    common = sorted(set(g570).intersection(g6244))
    out = config.COMMON_GENES
    out.write_text("\n".join(common) + "\n")
    print(f"GPL570 genes:                 {len(g570):>6}")
    print(f"GPL6244 genes:                {len(g6244):>6}")
    print(f"Common (LASSO candidate pool): {len(common):>6}")
    print(f"Wrote {out}")

    # Verify every gene in the expected panel is in the common pool
    missing = [g for g in config.EXPECTED_PANEL_GENES if g not in common]
    if missing:
        print(f"\nWARNING: expected panel genes missing from common pool: {missing}")
        print("(LASSO will not be able to select these as panel members.)")
    else:
        print(f"\nAll {len(config.EXPECTED_PANEL_GENES)} expected panel genes are in the common pool.")


if __name__ == "__main__":
    sys.exit(main() or 0)

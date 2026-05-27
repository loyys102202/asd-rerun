"""
11b_precompute_fold_change.py — Precompute bootstrap log2 fold change with 95 %
CI for each gene in the locked panel, on the training cohort GSE18123-GPL570.

This artifact is consumed by figures/20_make_sfig1.py (Panel C).

Inputs:
  - expr_gene_GPL570.parquet
  - phenotype_GPL570.csv
  - locked_model_main_7gene.json

Output:
  - panel_fold_change.json
"""
import json
import numpy as np
import pandas as pd

np.random.seed(42)
B = 2000  # bootstrap iterations

# ---- Load ----
expr = pd.read_parquet('expr_gene_GPL570.parquet')
pheno = pd.read_csv('phenotype_GPL570.csv')
with open('locked_model_main_7gene.json') as f:
    model = json.load(f)

# Handle gene × sample orientation
if not str(expr.index[0]).startswith('GSM'):
    expr = expr.T
sids = pheno['sample_id'].values
expr = expr.loc[sids]
y = pheno.set_index('sample_id').loc[expr.index, 'label'].astype(int).values

panel = model['panel']
results = {}
for g in panel:
    x = expr[g].astype(float).values
    asd = x[y == 1]
    ctl = x[y == 0]
    delta = asd.mean() - ctl.mean()        # log2 fold change on log2 scale
    deltas_b = np.empty(B)
    for b in range(B):
        a_b = np.random.choice(asd, size=len(asd), replace=True)
        c_b = np.random.choice(ctl, size=len(ctl), replace=True)
        deltas_b[b] = a_b.mean() - c_b.mean()
    results[g] = {
        'log2FC': float(delta),
        'CI_lo': float(np.percentile(deltas_b, 2.5)),
        'CI_hi': float(np.percentile(deltas_b, 97.5)),
        'FC_linear': float(2 ** delta),
        'mean_asd': float(asd.mean()),
        'mean_ctl': float(ctl.mean()),
        'n_asd': int(len(asd)),
        'n_ctl': int(len(ctl)),
    }
    print(f"  {g:10s}  log2FC = {delta:+.3f}  "
          f"[{results[g]['CI_lo']:+.3f}, {results[g]['CI_hi']:+.3f}]  "
          f"→ FC = {2**delta:.2f}")

with open('panel_fold_change.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Wrote panel_fold_change.json")

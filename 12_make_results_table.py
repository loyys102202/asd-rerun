"""12_make_results_table.py
================================================================================
Aggregate every numerical result from the pipeline into a single
results_summary.json (machine-readable) and Table_3.csv (paper-ready) for
direct cross-check against the Manuscript and Tables_revised.docx.

Reads everything in data/artifacts/ and produces a consolidated summary.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import config
from pipeline import read_json


def _safe_read(p: Path):
    if p.exists():
        try:
            return read_json(p)
        except json.JSONDecodeError:
            return None
    return None


def main():
    summary = {"manuscript_id": "ISCIENCE-D-26-03472",
                "panel_lock_threshold": config.STABILITY_THRESHOLD,
                "n_outer_folds_total": config.N_OUTER_FOLDS * config.N_REPEATS,
                "n_bootstrap_resamples": config.B_BOOTSTRAP,
                "n_permutations_training": config.N_PERMS_TRAIN,
                "n_permutations_external": config.N_PERMS_EXT}

    # --- locked panel ---
    if config.LOCKED_MODEL_MAIN.exists():
        m = read_json(config.LOCKED_MODEL_MAIN)
        summary["locked_panel"] = {
            "genes": m["genes"],
            "coef": dict(zip(m["genes"], m["coef"])),
            "intercept": m["intercept"],
            "train_mean": m["train_mean"],
            "train_sd": m["train_sd"],
            "apparent_auc": m.get("apparent_auc"),
            "apparent_brier": m.get("apparent_brier"),
        }
    if config.LOCKED_MODEL_REDUCED.exists():
        m = read_json(config.LOCKED_MODEL_REDUCED)
        summary["reduced_panel_for_GSE42133"] = {
            "genes": m["genes"],
            "coef": dict(zip(m["genes"], m["coef"])),
            "intercept": m["intercept"],
        }

    # --- selection frequencies ---
    if config.NESTEDCV_GENE_FREQ.exists():
        freq = pd.read_csv(config.NESTEDCV_GENE_FREQ)
        summary["selection_frequencies_top20"] = freq.head(20).to_dict(orient="records")

    # --- internal performance ---
    auc_by_rep = config.NESTEDCV_AUC_BY_REP
    if auc_by_rep.exists():
        rep_df = pd.read_csv(auc_by_rep)
        means = rep_df["mean_auc"].values
        summary["nested_cv"] = {
            "mean": float(means.mean()),
            "ci_2p5": float(np.percentile(means, 2.5)),
            "ci_97p5": float(np.percentile(means, 97.5)),
            "n_repeats": int(len(means)),
        }
    rc = _safe_read(config.ROBUSTNESS_CHECKS)
    if rc:
        summary["robustness_checks"] = rc

    # --- external validation: GPL6244 ---
    if config.EXTERNAL_GPL6244.exists():
        ext = pd.read_csv(config.EXTERNAL_GPL6244)
        summary["external_validation_GPL6244"] = ext.to_dict(orient="records")

    # --- external validation: GSE42133 ---
    summary["external_validation_GSE42133"] = _safe_read(config.EXTERNAL_GSE42133)

    # --- GSE25507 exploration ---
    summary["GSE25507_exploration"] = _safe_read(config.GSE25507_EXPLORATION)

    # --- calibration ---
    summary["calibration_GPL6244"] = _safe_read(config.CALIBRATION)

    # --- cell-marker correlations ---
    summary["cell_marker_correlations"] = _safe_read(config.CELL_MARKER_ANALYSIS)

    # Write
    from pipeline import write_json
    write_json(summary, config.RESULTS_SUMMARY)
    print(f"wrote {config.RESULTS_SUMMARY}")

    # Build a Table 3-equivalent CSV for paper cross-check
    rows = []

    # Internal
    if "locked_panel" in summary and summary["locked_panel"].get("apparent_auc") is not None:
        rows.append({"cohort": "GSE18123-GPL570", "tier": "Internal — apparent",
                     "n_ASD": 66, "n_ctrl": 33,
                     "auc": summary["locked_panel"]["apparent_auc"], "ci": "",
                     "brier": summary["locked_panel"]["apparent_brier"]})
    if "nested_cv" in summary:
        nc = summary["nested_cv"]
        rows.append({"cohort": "GSE18123-GPL570", "tier": "Internal — nested CV",
                     "n_ASD": 66, "n_ctrl": 33,
                     "auc": nc["mean"],
                     "ci": f"{nc['ci_2p5']:.3f}–{nc['ci_97p5']:.3f}",
                     "brier": ""})
    if "robustness_checks" in summary and "bootstrap_optimism" in summary["robustness_checks"]:
        bs = summary["robustness_checks"]["bootstrap_optimism"]
        rows.append({"cohort": "GSE18123-GPL570", "tier": "Internal — bootstrap (.632)",
                     "n_ASD": 66, "n_ctrl": 33,
                     "auc": bs["optimism_corrected_auc"], "ci": "",
                     "brier": ""})
    if "robustness_checks" in summary and "permutation_training" in summary["robustness_checks"]:
        pt = summary["robustness_checks"]["permutation_training"]
        rows.append({"cohort": "GSE18123-GPL570", "tier": "Internal — permutation null",
                     "n_ASD": 66, "n_ctrl": 33,
                     "auc": f"p = {pt['p_value']:.4f}", "ci": "", "brier": ""})

    # External GPL6244
    if "external_validation_GPL6244" in summary:
        for r in summary["external_validation_GPL6244"]:
            rows.append({
                "cohort": "GSE18123-GPL6244",
                "tier": r["analysis"],
                "n_ASD": r.get("n_ASD"),
                "n_ctrl": r.get("n_ctrl"),
                "auc": r.get("auc"),
                "ci": (f"{r.get('ci_low'):.3f}–{r.get('ci_high'):.3f}"
                       if pd.notna(r.get("ci_low")) and pd.notna(r.get("ci_high"))
                       else ""),
                "brier": r.get("brier"),
            })

    # External GSE42133
    if summary.get("external_validation_GSE42133"):
        gse42 = summary["external_validation_GSE42133"]
        for cfg_name, cfg in gse42.get("configurations", {}).items():
            rows.append({
                "cohort": "GSE42133",
                "tier": cfg_name,
                "n_ASD": gse42["n_ASD"], "n_ctrl": gse42["n_ctrl"],
                "auc": cfg["auc"],
                "ci": f"{cfg['ci_low']:.3f}–{cfg['ci_high']:.3f}",
                "brier": cfg["brier"],
            })

    # GSE25507
    if summary.get("GSE25507_exploration"):
        g25 = summary["GSE25507_exploration"]
        pa = g25.get("platform_adapted", {})
        if pa:
            rows.append({
                "cohort": "GSE25507", "tier": "Exploratory — platform-adapted",
                "n_ASD": g25["n_ASD"], "n_ctrl": g25["n_ctrl"],
                "auc": pa["auc"],
                "ci": f"{pa['ci_low']:.3f}–{pa['ci_high']:.3f}",
                "brier": "",
            })
        st = g25.get("strict_train_scale", {})
        if st:
            rows.append({
                "cohort": "GSE25507", "tier": "Exploratory — strict scale (sensitivity)",
                "n_ASD": g25["n_ASD"], "n_ctrl": g25["n_ctrl"],
                "auc": st["auc"], "ci": "", "brier": "",
            })
        wt = g25.get("within_tissue_refit", {})
        if wt:
            rows.append({
                "cohort": "GSE25507", "tier": "Exploratory — within-tissue refit",
                "n_ASD": g25["n_ASD"], "n_ctrl": g25["n_ctrl"],
                "auc": wt["refit_apparent_auc"], "ci": "", "brier": "",
            })

    df = pd.DataFrame(rows)
    out_csv = config.ARTIFACTS / "Table_3_reproduced.csv"
    df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")

    print("\nstage 12 complete.")
    print("\n--- summary ---")
    print(f"locked panel: {summary.get('locked_panel', {}).get('genes')}")
    nc = summary.get("nested_cv")
    if nc:
        print(f"nested CV AUC: {nc['mean']:.4f}  [{nc['ci_2p5']:.4f}, {nc['ci_97p5']:.4f}]")
    rc = summary.get("robustness_checks", {})
    if "bootstrap_optimism" in rc:
        bo = rc["bootstrap_optimism"]
        print(f"bootstrap-corrected AUC: {bo['optimism_corrected_auc']:.4f}  "
              f"(optimism {bo['mean_optimism']:.4f})")
    if "permutation_training" in rc:
        pt = rc["permutation_training"]
        print(f"training permutation p: {pt['p_value']:.4f}")


if __name__ == "__main__":
    sys.exit(main() or 0)

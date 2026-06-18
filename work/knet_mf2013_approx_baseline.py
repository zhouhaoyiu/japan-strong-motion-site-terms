#!/usr/bin/env python3
"""Approximate Morikawa-Fujiwara 2013 baseline for local K-NET targets.

This is a conservative bridge from the official J-SHIS/NIED MF2013 coefficient
CSV to the local 1997-2006 K-NET feature table. It is not a full MF2013
replication because the local K-NET cache lacks finite-fault distance,
mechanism/source-class labels, and moment magnitude for all events.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, median_absolute_error, r2_score


DEFAULT_WAVEFORM_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_record_waveform_features.csv"
)
DEFAULT_SA_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_record_response_spectrum_features.csv"
)
DEFAULT_SITE_MATCH = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/jshis_knet_site_schema_match.csv"
)
DEFAULT_COEFFS = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/work/external_data/jshis_mf2013/MF13rev_coefs.csv"
)
DEFAULT_EXISTING_METRICS = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_gmm_style_baseline_metrics.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)

TARGETS = [
    {
        "target": "log10_pga_horizontal_from_waveform_gal",
        "target_label": "Horizontal PGA",
        "period": "Acc.",
    },
    {
        "target": "log10_psa_0.3s_5pct_horizontal",
        "target_label": "SA(0.3s)",
        "period": "0.30",
    },
    {
        "target": "log10_psa_1s_5pct_horizontal",
        "target_label": "SA(1.0s)",
        "period": "1.00",
    },
    {
        "target": "log10_psa_3s_5pct_horizontal",
        "target_label": "SA(3.0s)",
        "period": "3.00",
    },
]

EVENT_TYPE_RULES = {
    "k1_all_crustal": 1,
    "k2_all_interface": 2,
    "k3_all_intraslab": 3,
    "depth_rule_k1_le30_else_k3": None,
}


def load_frame(waveform_path: Path, sa_path: Path, site_path: Path) -> pd.DataFrame:
    wave_cols = [
        "_id",
        "EventName",
        "Distance",
        "Magnitude",
        "Depth_km",
        "StationCode",
        "event_split",
        "time_split",
        "station_split",
        "log10_pga_horizontal_from_waveform_gal",
    ]
    wave = pd.read_csv(waveform_path, usecols=wave_cols, dtype={"_id": str, "EventName": str})
    sa_cols = [
        "_id",
        "log10_psa_0.3s_5pct_horizontal",
        "log10_psa_1s_5pct_horizontal",
        "log10_psa_3s_5pct_horizontal",
    ]
    sa = pd.read_csv(sa_path, usecols=sa_cols, dtype={"_id": str})
    site_cols = [
        "station_code",
        "jshis_match_found",
        "jshis_match_status",
        "jshis_avs30",
        "jshis_d1400",
        "jshis_dist_vf_mf13_nejapan",
        "jshis_dist_vf_mf13_swjapan",
    ]
    site = pd.read_csv(site_path, usecols=site_cols)
    site = site[site["jshis_match_found"].eq(True)].drop_duplicates("station_code")
    df = wave.merge(sa, on="_id", how="left").merge(
        site,
        left_on="StationCode",
        right_on="station_code",
        how="inner",
    )
    for col in [
        "Distance",
        "Magnitude",
        "Depth_km",
        "jshis_avs30",
        "jshis_d1400",
        "jshis_dist_vf_mf13_nejapan",
        "jshis_dist_vf_mf13_swjapan",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[
        df["Distance"].notna()
        & df["Magnitude"].notna()
        & df["Depth_km"].notna()
        & df["jshis_avs30"].notna()
    ].copy()
    return df


def load_coeffs(path: Path) -> pd.DataFrame:
    coeffs = pd.read_csv(path)
    coeffs.columns = [col.strip() for col in coeffs.columns]
    coeffs["Period"] = coeffs["Period"].astype(str).str.strip()
    for col in coeffs.columns:
        if col != "Period":
            coeffs[col] = pd.to_numeric(coeffs[col], errors="coerce")
    return coeffs.set_index("Period")


def source_class(rule_name: str, depth_km: pd.Series) -> pd.Series:
    fixed = EVENT_TYPE_RULES[rule_name]
    if fixed is not None:
        return pd.Series(fixed, index=depth_km.index)
    return pd.Series(np.where(depth_km <= 30.0, 1, 3), index=depth_km.index)


def mf2013_prediction(
    df: pd.DataFrame,
    coeff: pd.Series,
    rule_name: str,
    site_terms: bool,
) -> pd.Series:
    mw_proxy = np.minimum(df["Magnitude"].astype(float), 8.2)
    x = np.maximum(df["Distance"].astype(float), 1.0)
    k = source_class(rule_name, df["Depth_km"].astype(float))
    b = np.select([k.eq(1), k.eq(2), k.eq(3)], [coeff["b1"], coeff["b2"], coeff["b3"]], default=coeff["b1"])
    c = np.select([k.eq(1), k.eq(2), k.eq(3)], [coeff["c1"], coeff["c2"], coeff["c3"]], default=coeff["c1"])
    log10_a = (
        coeff["a"] * (mw_proxy - 16.0) ** 2
        + b * x
        + c
        - np.log10(x + coeff["d"] * 10.0 ** (0.5 * mw_proxy))
    )
    if site_terms:
        d1400 = df["jshis_d1400"].astype(float)
        avs30 = df["jshis_avs30"].astype(float)
        gd = np.where(
            d1400.notna() & (d1400 > 0),
            coeff["pd"] * np.log10(np.maximum(coeff["Dlmin"], d1400) / 300.0),
            0.0,
        )
        gs = np.where(
            avs30.notna() & (avs30 > 0),
            coeff["ps"] * np.log10(np.minimum(coeff["Vsmax"], avs30) / 350.0),
            0.0,
        )
        log10_a = log10_a + gd + gs
    return pd.Series(log10_a, index=df.index)


def split_masks(df: pd.DataFrame, split_col: str) -> tuple[pd.Series, pd.Series]:
    train = df[split_col].isin(["train", "dev"])
    test = df[split_col].eq("test")
    return train, test


def calibrate(train_pred: np.ndarray, train_y: np.ndarray, test_pred: np.ndarray) -> tuple[np.ndarray, float, float]:
    mask = np.isfinite(train_pred) & np.isfinite(train_y)
    if mask.sum() < 10 or np.nanstd(train_pred[mask]) < 1e-9:
        intercept = float(np.nanmean(train_y[mask])) if mask.sum() else 0.0
        slope = 0.0
    else:
        slope, intercept = np.polyfit(train_pred[mask], train_y[mask], deg=1)
        slope = float(slope)
        intercept = float(intercept)
    return intercept + slope * test_pred, intercept, slope


def metric_row(
    split: str,
    target: str,
    target_label: str,
    model: str,
    train_rows: int,
    test_rows: int,
    y_test: np.ndarray,
    pred_test: np.ndarray,
    calibration_intercept: float | None,
    calibration_slope: float | None,
) -> dict:
    mask = np.isfinite(y_test) & np.isfinite(pred_test)
    y = y_test[mask]
    pred = pred_test[mask]
    residual = y - pred
    return {
        "split": split,
        "target": target,
        "target_label": target_label,
        "model": model,
        "train_rows": train_rows,
        "test_rows": test_rows,
        "usable_test_rows": int(mask.sum()),
        "mae": float(mean_absolute_error(y, pred)),
        "median_ae": float(median_absolute_error(y, pred)),
        "r2": float(r2_score(y, pred)),
        "bias_observed_minus_pred": float(np.mean(residual)),
        "calibration_intercept": calibration_intercept,
        "calibration_slope": calibration_slope,
    }


def evaluate(df: pd.DataFrame, coeffs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    residual_frames = []
    for target_info in TARGETS:
        target = target_info["target"]
        label = target_info["target_label"]
        coeff = coeffs.loc[target_info["period"]]
        work = df[df[target].notna()].copy()
        y_all = work[target].astype(float).to_numpy()
        for rule_name in EVENT_TYPE_RULES:
            for site_terms in [False, True]:
                pred_all = mf2013_prediction(work, coeff, rule_name, site_terms).to_numpy()
                site_label = "site" if site_terms else "basic"
                base_name = f"mf2013_{site_label}_{rule_name}"
                for split in ["event_split", "time_split", "station_split"]:
                    train_mask, test_mask = split_masks(work, split)
                    train_idx = np.where(train_mask.to_numpy())[0]
                    test_idx = np.where(test_mask.to_numpy())[0]
                    if len(train_idx) == 0 or len(test_idx) == 0:
                        continue
                    y_test = y_all[test_idx]
                    pred_test_raw = pred_all[test_idx]
                    metric_rows.append(
                        metric_row(
                            split,
                            target,
                            label,
                            base_name + "_raw",
                            int(len(train_idx)),
                            int(len(test_idx)),
                            y_test,
                            pred_test_raw,
                            None,
                            None,
                        )
                    )

                    pred_test_cal, intercept, slope = calibrate(
                        pred_all[train_idx],
                        y_all[train_idx],
                        pred_test_raw,
                    )
                    metric_rows.append(
                        metric_row(
                            split,
                            target,
                            label,
                            base_name + "_train_calibrated",
                            int(len(train_idx)),
                            int(len(test_idx)),
                            y_test,
                            pred_test_cal,
                            intercept,
                            slope,
                        )
                    )
                    residual_frames.append(
                        pd.DataFrame(
                            {
                                "EventName": work.iloc[test_idx]["EventName"].to_numpy(),
                                "StationCode": work.iloc[test_idx]["StationCode"].to_numpy(),
                                "split": split,
                                "target": target,
                                "target_label": label,
                                "model": base_name + "_train_calibrated",
                                "observed_log10": y_test,
                                "predicted_log10": pred_test_cal,
                                "residual_log10": y_test - pred_test_cal,
                            }
                        )
                    )
    return pd.DataFrame(metric_rows), pd.concat(residual_frames, ignore_index=True)


def compare_existing(mf_metrics: pd.DataFrame, existing_metrics_path: Path) -> pd.DataFrame:
    existing = pd.read_csv(existing_metrics_path)
    existing = existing[existing["target_label"].isin([item["target_label"] for item in TARGETS])]
    rows = []
    for (split, target_label), sub in existing.groupby(["split", "target_label"]):
        mf_sub = mf_metrics[(mf_metrics["split"] == split) & (mf_metrics["target_label"] == target_label)]
        if mf_sub.empty:
            continue
        raw = mf_sub[mf_sub["model"].str.endswith("_raw")].sort_values("mae").iloc[0]
        cal = mf_sub[mf_sub["model"].str.endswith("_train_calibrated")].sort_values("mae").iloc[0]
        gmm = sub[sub["model"].eq("gmm_style")].iloc[0]
        best = sub.sort_values("mae").iloc[0]
        rows.append(
            {
                "split": split,
                "target_label": target_label,
                "best_mf2013_raw_model": raw["model"],
                "best_mf2013_raw_mae": raw["mae"],
                "best_mf2013_calibrated_model": cal["model"],
                "best_mf2013_calibrated_mae": cal["mae"],
                "empirical_gmm_style_mae": gmm["mae"],
                "best_existing_model": best["model"],
                "best_existing_mae": best["mae"],
                "best_existing_reduction_vs_calibrated_mf2013_pct": (
                    (cal["mae"] - best["mae"]) / cal["mae"] * 100.0
                ),
                "gmm_style_reduction_vs_calibrated_mf2013_pct": (
                    (cal["mae"] - gmm["mae"]) / cal["mae"] * 100.0
                ),
            }
        )
    return pd.DataFrame(rows)


def filter_best_calibrated_residuals(residuals: pd.DataFrame, comparison: pd.DataFrame) -> pd.DataFrame:
    keep = comparison[["split", "target_label", "best_mf2013_calibrated_model"]].rename(
        columns={"best_mf2013_calibrated_model": "model"}
    )
    return residuals.merge(keep, on=["split", "target_label", "model"], how="inner")


def plot_comparison(comparison: pd.DataFrame, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    order = ["Horizontal PGA", "SA(0.3s)", "SA(1.0s)", "SA(3.0s)"]
    for split, sub in comparison.groupby("split"):
        sub = sub.set_index("target_label").loc[[label for label in order if label in sub["target_label"].values]]
        values = sub[
            [
                "best_mf2013_calibrated_mae",
                "empirical_gmm_style_mae",
                "best_existing_mae",
            ]
        ]
        x = np.arange(len(values.index))
        width = 0.25
        fig, ax = plt.subplots(figsize=(8.2, 4.6))
        ax.bar(x - width, values["best_mf2013_calibrated_mae"], width, label="MF2013 approx calibrated")
        ax.bar(x, values["empirical_gmm_style_mae"], width, label="Empirical GMM-style")
        ax.bar(x + width, values["best_existing_mae"], width, label="Best existing")
        ax.set_xticks(x)
        ax.set_xticklabels(values.index, rotation=15, ha="right")
        ax.set_ylabel("MAE in log10 units")
        ax.set_title(f"MF2013 approximate baseline comparison: {split}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(figure_dir / f"knet_mf2013_approx_baseline_comparison_{split}.png", dpi=180)
        plt.close(fig)


def write_report(metrics: pd.DataFrame, comparison: pd.DataFrame, output_dir: Path, records_total: int) -> None:
    lines = ["# K-NET Approximate Morikawa-Fujiwara 2013 Baseline\n\n"]
    lines.append("## Source and implementation\n")
    lines.append(
        "- Coefficients: official J-SHIS/NIED `MF13rev_coefs.csv`, downloaded from the Morikawa and Fujiwara (2013) GMPE page.\n"
    )
    lines.append("- Source URL: `https://www.j-shis.bosai.go.jp/en/labs/mf2013/`.\n")
    lines.append(
        "- Formula used: MF2013 basic magnitude-distance term plus optional J-SHIS `D1400` and `AVS30` site corrections.\n"
    )
    lines.append("- Correction terms not used in this local approximate baseline: anomalous intensity `AI` and Philippine Sea Plate `PH`.\n")
    lines.append("- Local substitutions: K-NET `Magnitude` is used as an Mw proxy, and K-NET `Distance` is used as the closest available distance proxy.\n")
    lines.append("- Source-class uncertainty is handled by evaluating fixed k=1/k=2/k=3 and a depth rule k=1 for depth <=30 km else k=3.\n\n")

    lines.append("## Usable data\n")
    lines.append(f"- Records after official site matching and required input filtering: {records_total:,}.\n")
    lines.append("- Targets evaluated: horizontal waveform PGA and 5%-damped SA at 0.3, 1.0, and 3.0 s.\n\n")

    lines.append("## Best calibrated MF2013 approximate baseline vs existing models\n")
    for split in ["event_split", "time_split", "station_split"]:
        lines.append(f"### {split}\n")
        sub = comparison[comparison["split"] == split]
        for _, row in sub.iterrows():
            lines.append(
                f"- {row['target_label']}: best MF2013 calibrated MAE={row['best_mf2013_calibrated_mae']:.4f} "
                f"({row['best_mf2013_calibrated_model']}), empirical GMM-style MAE={row['empirical_gmm_style_mae']:.4f}, "
                f"best existing={row['best_existing_model']} MAE={row['best_existing_mae']:.4f}, "
                f"best-existing reduction vs calibrated MF2013={row['best_existing_reduction_vs_calibrated_mf2013_pct']:.1f}%.\n"
            )
        lines.append("\n")

    lines.append("## Raw MF2013 check\n")
    raw_best = (
        metrics[metrics["model"].str.endswith("_raw")]
        .sort_values(["split", "target_label", "mae"])
        .groupby(["split", "target_label"])
        .head(1)
    )
    for _, row in raw_best.iterrows():
        lines.append(
            f"- {row['split']} / {row['target_label']}: best raw={row['model']} "
            f"MAE={row['mae']:.4f}, bias observed-minus-pred={row['bias_observed_minus_pred']:.4f}.\n"
        )

    lines.append("\n## Outputs\n")
    for name in [
        "knet_mf2013_approx_baseline_metrics.csv",
        "knet_mf2013_approx_baseline_comparison.csv",
        "knet_mf2013_approx_baseline_best_residuals.csv",
        "figures/knet_mf2013_approx_baseline_comparison_event_split.png",
        "figures/knet_mf2013_approx_baseline_comparison_time_split.png",
        "figures/knet_mf2013_approx_baseline_comparison_station_split.png",
    ]:
        lines.append(f"- `{name}`\n")

    lines.append("\n## Guardrails\n")
    lines.append(
        "- This is a published-GMM approximate benchmark, not a complete MF2013 reproduction. It should not be presented as exact GMPE residuals.\n"
    )
    lines.append(
        "- The calibrated scores are valid for model comparison but are no longer pure published-GMM predictions because the calibration uses the local train split.\n"
    )
    lines.append(
        "- A manuscript should either implement the full MF2013 input requirements or use the official J-SHIS/NIED `smrec_schema.tsv` records for a separate official-data residual analysis.\n"
    )
    (output_dir / "knet_mf2013_approx_baseline.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--waveform-input", type=Path, default=DEFAULT_WAVEFORM_INPUT)
    parser.add_argument("--sa-input", type=Path, default=DEFAULT_SA_INPUT)
    parser.add_argument("--site-match", type=Path, default=DEFAULT_SITE_MATCH)
    parser.add_argument("--coeffs", type=Path, default=DEFAULT_COEFFS)
    parser.add_argument("--existing-metrics", type=Path, default=DEFAULT_EXISTING_METRICS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_frame(args.waveform_input, args.sa_input, args.site_match)
    coeffs = load_coeffs(args.coeffs)
    metrics, residuals = evaluate(df, coeffs)
    comparison = compare_existing(metrics, args.existing_metrics)
    best_residuals = filter_best_calibrated_residuals(residuals, comparison)

    metrics.to_csv(args.output_dir / "knet_mf2013_approx_baseline_metrics.csv", index=False)
    comparison.to_csv(args.output_dir / "knet_mf2013_approx_baseline_comparison.csv", index=False)
    best_residuals.to_csv(args.output_dir / "knet_mf2013_approx_baseline_best_residuals.csv", index=False)
    old_residuals = args.output_dir / "knet_mf2013_approx_baseline_residuals.csv"
    if old_residuals.exists():
        old_residuals.unlink()
    plot_comparison(comparison, args.output_dir)
    write_report(metrics, comparison, args.output_dir, records_total=len(df))

    print(f"records_after_site_match={len(df):,}")
    print(f"metrics_rows={len(metrics):,}")
    print(f"comparison_rows={len(comparison):,}")
    print(f"best_residual_rows={len(best_residuals):,}")


if __name__ == "__main__":
    main()

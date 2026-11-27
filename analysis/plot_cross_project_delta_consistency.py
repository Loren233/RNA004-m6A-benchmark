
import os as _os
from pathlib import Path as _Path
RNA004_STAT_ROOT = _os.environ.get(
    "RNA004_STAT_ROOT",
    str(_Path(__file__).resolve().parents[0])
)

#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


ROOT = Path(f"{RNA004_STAT_ROOT}")

DIR_106788 = ROOT / "106788" / "delta_beta"
DIR_80229 = ROOT / "80229" / "delta_beta"

TOOLS = ["deeprm", "dorado", "m6anet", "singlemod"]

OUTDIR = ROOT / "cross_project_delta_consistency"
OUTDIR.mkdir(parents=True, exist_ok=True)

MERGED_DIR = OUTDIR / "merged_tables"
MERGED_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_ROWS = []

PLOT_SAMPLE_MAX = 200000


def load_delta_table(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = ["chr", "pos_1based", "strand", "delta_beta"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"{path} missing columns: {missing}")

    df["chr"] = df["chr"].astype("string").str.strip()
    df["strand"] = df["strand"].astype("string").str.strip()
    df["pos_1based"] = pd.to_numeric(df["pos_1based"], errors="coerce")
    df["delta_beta"] = pd.to_numeric(df["delta_beta"], errors="coerce")

    # 可选列
    if "min_depth" in df.columns:
        df["min_depth"] = pd.to_numeric(df["min_depth"], errors="coerce")
    else:
        df["min_depth"] = np.nan

    df = df.dropna(subset=["chr", "pos_1based", "strand", "delta_beta"]).copy()
    df["pos_1based"] = df["pos_1based"].astype(int)

    # 只保留合法delta_beta
    df = df[
        (df["delta_beta"] >= -1) &
        (df["delta_beta"] <= 1)
    ].copy()

    df["site_id"] = (
        df["chr"].astype(str) + "|" +
        df["pos_1based"].astype(str) + "|" +
        df["strand"].astype(str)
    )

    # 保留需要的列
    df = df[
        ["chr", "pos_1based", "strand", "site_id", "delta_beta", "min_depth"]
    ].copy()

    # 防御性检查：delta_beta文件本身不应重复
    dup = df.duplicated(["chr", "pos_1based", "strand"], keep=False)
    if dup.any():
        ndup = int(dup.sum())
        raise RuntimeError(f"{path} still contains {ndup} duplicated sites")

    return df


def compute_sign_concordance(x, y):
    sign_x = np.sign(x)
    sign_y = np.sign(y)
    return np.mean(sign_x == sign_y)


def make_plot(df, tool, outpath):
    n = len(df)

    x = df["delta_beta_106788"].to_numpy()
    y = df["delta_beta_80229"].to_numpy()

    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)

    mae = np.mean(np.abs(y - x))
    mean_diff = np.mean(y - x)
    sign_concordance = compute_sign_concordance(x, y)

    quadrant_same = np.mean(
        ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
    )

    # 作图抽样
    if n > PLOT_SAMPLE_MAX:
        plot_df = df.sample(PLOT_SAMPLE_MAX, random_state=42)
    else:
        plot_df = df

    x_plot = plot_df["delta_beta_106788"].to_numpy()
    y_plot = plot_df["delta_beta_80229"].to_numpy()

    plt.figure(figsize=(6.8, 6.8))

    plt.scatter(
        x_plot,
        y_plot,
        s=4,
        alpha=0.20,
        rasterized=True
    )

    # 参考线
    plt.plot([-1, 1], [-1, 1], linestyle="--", linewidth=1)
    plt.axhline(0, linestyle=":", linewidth=1)
    plt.axvline(0, linestyle=":", linewidth=1)

    plt.xlim(-1, 1)
    plt.ylim(-1, 1)

    plt.xlabel("Δβ_106788 = β_HEK293T - β_HeLa")
    plt.ylabel("Δβ_80229 = β_HEK293T - β_HeLa")
    plt.title(f"{tool} | Δβ_80229 vs Δβ_106788")

    text = (
        f"n = {n:,}\n"
        f"Pearson r = {pearson_r:.3f}\n"
        f"Spearman ρ = {spearman_rho:.3f}\n"
        f"MAE = {mae:.3f}\n"
        f"Mean(80229 - 106788) = {mean_diff:.3f}\n"
        f"Sign concordance = {sign_concordance:.3f}\n"
        f"Same-direction quadrants = {quadrant_same:.3f}"
    )

    plt.text(
        0.03, 0.97, text,
        transform=plt.gca().transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", alpha=0.15)
    )

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

    return {
        "tool": tool,
        "n_sites": n,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
        "mae_delta_beta_80229_vs_106788": mae,
        "mean_delta_beta_80229_minus_106788": mean_diff,
        "sign_concordance_raw": sign_concordance,
        "same_direction_quadrant_fraction": quadrant_same,
        "plot_file": str(outpath),
    }


def main():
    for tool in TOOLS:
        file_106788 = DIR_106788 / f"{tool}.HEK_minus_HeLa.delta_beta.tsv"
        file_80229 = DIR_80229 / f"{tool}.HEK_minus_HeLa.delta_beta.tsv"

        if not file_106788.exists():
            print(f"[WARNING] missing file: {file_106788}")
            continue
        if not file_80229.exists():
            print(f"[WARNING] missing file: {file_80229}")
            continue

        print(f"[INFO] processing tool: {tool}")

        df_106788 = load_delta_table(file_106788)
        df_80229 = load_delta_table(file_80229)

        df_106788 = df_106788.rename(columns={
            "delta_beta": "delta_beta_106788",
            "min_depth": "min_depth_106788",
        })

        df_80229 = df_80229.rename(columns={
            "delta_beta": "delta_beta_80229",
            "min_depth": "min_depth_80229",
        })

        merged = df_106788.merge(
            df_80229,
            on=["chr", "pos_1based", "strand", "site_id"],
            how="inner",
            validate="one_to_one"
        )

        merged["abs_diff_delta_beta"] = np.abs(
            merged["delta_beta_80229"] - merged["delta_beta_106788"]
        )

        merged["direction_106788"] = np.select(
            [
                merged["delta_beta_106788"] > 0,
                merged["delta_beta_106788"] < 0,
            ],
            [
                "HEK_higher",
                "HeLa_higher",
            ],
            default="equal",
        )

        merged["direction_80229"] = np.select(
            [
                merged["delta_beta_80229"] > 0,
                merged["delta_beta_80229"] < 0,
            ],
            [
                "HEK_higher",
                "HeLa_higher",
            ],
            default="equal",
        )

        merged["direction_concordant_raw"] = (
            np.sign(merged["delta_beta_106788"]) ==
            np.sign(merged["delta_beta_80229"])
        )

        merged_out = MERGED_DIR / f"{tool}.cross_project.delta_beta.merged.tsv"
        merged.to_csv(merged_out, sep="\t", index=False)

        plot_out = OUTDIR / f"{tool}.delta_beta_80229_vs_106788.png"
        row = make_plot(merged, tool, plot_out)
        row["merged_table"] = str(merged_out)

        SUMMARY_ROWS.append(row)

    summary_df = pd.DataFrame(SUMMARY_ROWS)
    summary_path = OUTDIR / "cross_project_delta_consistency_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print(f"[INFO] summary written to: {summary_path}")
    print(f"[INFO] plots written to: {OUTDIR}")


if __name__ == "__main__":
    main()


import os as _os
from pathlib import Path as _Path
RNA004_STAT_ROOT = _os.environ.get(
    "RNA004_STAT_ROOT",
    str(_Path(__file__).resolve().parents[0])
)

#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


ROOT = Path(f"{RNA004_STAT_ROOT}")

PROJECTS = {
    "106788": ROOT / "106788" / "delta_beta",
    "80229": ROOT / "80229" / "delta_beta",
}

TOOLS = ["deeprm", "dorado", "m6anet", "singlemod"]

OUTDIR = ROOT / "beta_consistency_plots"
OUTDIR.mkdir(parents=True, exist_ok=True)

SUMMARY_ROWS = []

# 为了避免上千万点直接scatter太卡，这里只抽样用于作图
PLOT_SAMPLE_MAX = 200000


def load_delta_table(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = ["beta_hek", "beta_hela"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"{path} missing columns: {missing}")

    df["beta_hek"] = pd.to_numeric(df["beta_hek"], errors="coerce")
    df["beta_hela"] = pd.to_numeric(df["beta_hela"], errors="coerce")

    df = df.dropna(subset=["beta_hek", "beta_hela"]).copy()

    # 只保留合法beta
    df = df[
        (df["beta_hek"] >= 0) & (df["beta_hek"] <= 1) &
        (df["beta_hela"] >= 0) & (df["beta_hela"] <= 1)
    ].copy()

    return df


def make_plot(df, project, tool, outpath):
    n = len(df)

    x = df["beta_hek"].to_numpy()
    y = df["beta_hela"].to_numpy()

    # 统计量用全数据
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)

    mean_abs_diff = np.mean(np.abs(x - y))
    mean_diff = np.mean(x - y)

    # 画图只抽样，避免太慢
    if n > PLOT_SAMPLE_MAX:
        plot_df = df.sample(PLOT_SAMPLE_MAX, random_state=42)
    else:
        plot_df = df

    x_plot = plot_df["beta_hek"].to_numpy()
    y_plot = plot_df["beta_hela"].to_numpy()

    plt.figure(figsize=(6.5, 6.5))

    plt.scatter(
        x_plot,
        y_plot,
        s=4,
        alpha=0.20,
        rasterized=True
    )

    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)

    plt.xlim(0, 1)
    plt.ylim(0, 1)

    plt.xlabel("β_HEK293T")
    plt.ylabel("β_HeLa")
    plt.title(f"{project} | {tool} | β_HEK293T vs β_HeLa")

    text = (
        f"n = {n:,}\n"
        f"Pearson r = {pearson_r:.3f}\n"
        f"Spearman ρ = {spearman_rho:.3f}\n"
        f"Mean(β_HEK - β_HeLa) = {mean_diff:.3f}\n"
        f"MAE = {mean_abs_diff:.3f}"
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
        "project": project,
        "tool": tool,
        "n_sites": n,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
        "mean_beta_hek_minus_hela": mean_diff,
        "mae_beta_hek_vs_hela": mean_abs_diff,
        "plot_file": str(outpath),
    }


def main():
    for project, project_dir in PROJECTS.items():
        for tool in TOOLS:
            infile = project_dir / f"{tool}.HEK_minus_HeLa.delta_beta.tsv"
            if not infile.exists():
                print(f"[WARNING] missing file: {infile}")
                continue

            print(f"[INFO] processing {infile}")
            df = load_delta_table(infile)

            outpath = OUTDIR / f"{project}.{tool}.beta_hek_vs_beta_hela.png"

            row = make_plot(df, project, tool, outpath)
            SUMMARY_ROWS.append(row)

    summary_df = pd.DataFrame(SUMMARY_ROWS)
    summary_path = OUTDIR / "beta_consistency_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print(f"[INFO] summary written to: {summary_path}")
    print(f"[INFO] plots written to: {OUTDIR}")


if __name__ == "__main__":
    main()

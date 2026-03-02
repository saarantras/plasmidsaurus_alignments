#!/usr/bin/env python3
"""Summarize per-feature coverage using BED annotations and samtools depth output."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, help="BED file with 6 columns")
    parser.add_argument("--depth", required=True, help="samtools depth -aa output")
    parser.add_argument("--output", required=True, help="Output TSV summary")
    parser.add_argument("--plot", help="Optional output barplot PNG")
    return parser.parse_args()


def load_features(path: Path) -> pd.DataFrame:
    cols = ["contig", "start", "end", "feature_type", "feature_id", "strand"]
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.empty:
        return pd.DataFrame(columns=cols)

    if df.shape[1] < 3:
        raise ValueError(f"Feature BED needs at least 3 columns: {path}")

    while df.shape[1] < 6:
        df[df.shape[1]] = "."

    df = df.iloc[:, :6]
    df.columns = cols
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def load_depth(path: Path) -> pd.DataFrame:
    cols = ["contig", "pos", "depth"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    if df.empty:
        return df
    df["pos"] = df["pos"].astype(int)
    df["depth"] = df["depth"].astype(float)
    return df


def summarize(features: pd.DataFrame, depth: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, feat in features.iterrows():
        contig = feat["contig"]
        start = int(feat["start"])
        end = int(feat["end"])
        length = max(0, end - start)

        # BED is 0-based, end-exclusive; samtools depth is 1-based.
        mask = (
            (depth["contig"] == contig)
            & (depth["pos"] >= start + 1)
            & (depth["pos"] <= end)
        )
        segment = depth.loc[mask, "depth"]

        if length == 0 or segment.empty:
            covered_bases = 0
            mean_depth = 0.0
            median_depth = 0.0
            max_depth = 0.0
            total_depth = 0.0
        else:
            covered_bases = int((segment > 0).sum())
            mean_depth = float(segment.mean())
            median_depth = float(segment.median())
            max_depth = float(segment.max())
            total_depth = float(segment.sum())

        frac_covered = (covered_bases / length) if length > 0 else 0.0

        rows.append(
            {
                "contig": contig,
                "start": start,
                "end": end,
                "feature_type": feat["feature_type"],
                "feature_id": feat["feature_id"],
                "strand": feat["strand"],
                "length": length,
                "covered_bases": covered_bases,
                "frac_bases_covered": frac_covered,
                "mean_depth": mean_depth,
                "median_depth": median_depth,
                "max_depth": max_depth,
                "total_depth": total_depth,
            }
        )

    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values(["start", "end", "feature_type", "feature_id"])
    return summary


def plot_summary(summary: pd.DataFrame, plot_path: Path) -> None:
    if summary.empty:
        return

    labels = [f"{t}:{i}" for t, i in zip(summary["feature_type"], summary["feature_id"])]
    x = range(len(summary))

    width = max(9, int(0.7 * len(summary)))
    plt.figure(figsize=(width, 4.8))
    plt.bar(x, summary["mean_depth"], color="#4C78A8")
    plt.xticks(list(x), labels, rotation=65, ha="right", fontsize=8)
    plt.ylabel("Mean depth")
    plt.title("Mean coverage depth by annotated feature")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    plt.close()


def main() -> None:
    args = parse_args()
    feature_path = Path(args.features)
    depth_path = Path(args.depth)
    output_path = Path(args.output)

    features = load_features(feature_path)
    depth = load_depth(depth_path)

    summary = summarize(features, depth)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote feature summary: {output_path} ({len(summary)} features)")

    if args.plot:
        plot_path = Path(args.plot)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        plot_summary(summary, plot_path)
        print(f"Wrote feature summary plot: {plot_path}")


if __name__ == "__main__":
    main()

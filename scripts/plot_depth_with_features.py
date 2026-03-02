#!/usr/bin/env python3
"""Plot per-base depth with annotation tracks from BED features."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", required=True, help="samtools depth -aa output")
    parser.add_argument("--features", required=True, help="Feature BED file")
    parser.add_argument("--output", required=True, help="Output PNG")
    parser.add_argument("--title", default="Depth with annotation tracks")
    return parser.parse_args()


def load_depth(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=["contig", "pos", "depth"])
    if df.empty:
        raise ValueError(f"Depth file is empty: {path}")
    df["pos"] = df["pos"].astype(int)
    df["depth"] = df["depth"].astype(float)
    return df


def load_features(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.empty:
        return pd.DataFrame(columns=["contig", "start", "end", "feature_type", "feature_id", "strand"])

    while df.shape[1] < 6:
        df[df.shape[1]] = "."
    df = df.iloc[:, :6]
    df.columns = ["contig", "start", "end", "feature_type", "feature_id", "strand"]
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def main() -> None:
    args = parse_args()
    depth = load_depth(Path(args.depth))
    features = load_features(Path(args.features))

    contig = depth["contig"].iloc[0]
    depth = depth[depth["contig"] == contig]
    features = features[features["contig"] == contig]

    feature_types = sorted(features["feature_type"].unique().tolist()) if not features.empty else []
    y_map = {ft: i for i, ft in enumerate(feature_types)}
    cmap = plt.get_cmap("tab20")
    color_map = {ft: cmap(i % 20) for i, ft in enumerate(feature_types)}

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(12, 6),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax1.plot(depth["pos"], depth["depth"], color="#1f4e79", linewidth=1.0)
    ax1.set_ylabel("Depth")
    ax1.set_title(args.title)
    ax1.grid(alpha=0.2)

    if features.empty:
        ax2.text(0.5, 0.5, "No features", transform=ax2.transAxes, ha="center", va="center")
        ax2.set_yticks([])
    else:
        for _, row in features.iterrows():
            x_start = row["start"] + 1
            width = max(1, row["end"] - row["start"])
            y = y_map[row["feature_type"]]
            ax2.broken_barh([(x_start, width)], (y - 0.35, 0.7), facecolors=color_map[row["feature_type"]])
            if width >= 18:
                ax2.text(
                    x_start + width / 2,
                    y,
                    str(row["feature_id"]),
                    ha="center",
                    va="center",
                    fontsize=7,
                    clip_on=True,
                )

        ax2.set_ylim(-0.8, len(feature_types) - 0.2)
        ax2.set_yticks(list(y_map.values()))
        ax2.set_yticklabels(list(y_map.keys()))

    ax2.set_xlabel(f"Position on {contig}")
    ax2.grid(alpha=0.2, axis="x")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Wrote depth+feature plot: {out_path}")


if __name__ == "__main__":
    main()

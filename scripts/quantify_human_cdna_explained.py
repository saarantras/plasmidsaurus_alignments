#!/usr/bin/env python3
"""Quantify how much of a read set is explainable as human RNA cDNA.

Transcriptome-only workflow:
1) align reads to human transcriptome with minimap2 (map-ont)
2) pick a best alignment per read
3) call reads "explained" if best query coverage and identity pass thresholds
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, help="FASTQ/FASTQ.GZ reads")
    parser.add_argument("--ref-mmi", required=True, help="minimap2 transcriptome index (.mmi)")
    parser.add_argument("--threads", type=int, default=4, help="minimap2 threads")
    parser.add_argument("--min-query-cov", type=float, default=0.50, help="Explained threshold")
    parser.add_argument("--min-identity", type=float, default=0.80, help="Explained threshold")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--keep-paf",
        action="store_true",
        help="Keep uncompressed PAF in addition to .paf.gz",
    )
    return parser.parse_args()


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def iter_fastq_lengths(path: Path) -> Iterator[Tuple[str, int]]:
    with open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            yield header[1:].split()[0], len(seq)


def parse_gencode_target_name(target: str) -> Dict[str, str]:
    parts = target.split("|")
    out = {
        "transcript_id": parts[0] if len(parts) > 0 else target,
        "gene_id": parts[1] if len(parts) > 1 else "",
        "transcript_name": parts[4] if len(parts) > 4 else "",
        "gene_name": parts[5] if len(parts) > 5 else "",
        "biotype": parts[7] if len(parts) > 7 else "",
    }
    return out


def ensure_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Missing required executable in PATH: {tool}")


def run_minimap2_paf(reads: Path, ref_mmi: Path, threads: int, paf_out: Path) -> None:
    cmd = [
        "minimap2",
        "-x",
        "map-ont",
        "-t",
        str(threads),
        str(ref_mmi),
        str(reads),
    ]
    with paf_out.open("w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)


def select_best_alignment(df: pd.DataFrame) -> pd.DataFrame:
    # Maximize query coverage first, then identity, then nmatch.
    df = df.copy()
    df["sort_key_1"] = df["query_cov"]
    df["sort_key_2"] = df["identity"]
    df["sort_key_3"] = df["nmatch"]
    df = df.sort_values(
        ["read_id", "sort_key_1", "sort_key_2", "sort_key_3"],
        ascending=[True, False, False, False],
    )
    best = df.groupby("read_id", as_index=False).first()
    return best.drop(columns=["sort_key_1", "sort_key_2", "sort_key_3"])


def build_threshold_grid(best_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cov_thresholds = [round(x / 100.0, 2) for x in range(30, 95, 5)]
    id_thresholds = [round(x / 100.0, 2) for x in range(70, 100, 5)]
    n = len(best_df)
    for cov_t in cov_thresholds:
        for id_t in id_thresholds:
            explained = (best_df["best_query_cov"] >= cov_t) & (best_df["best_identity"] >= id_t)
            cnt = int(explained.sum())
            rows.append(
                {
                    "query_cov_threshold": cov_t,
                    "identity_threshold": id_t,
                    "explained_reads": cnt,
                    "explained_fraction": cnt / n if n else 0.0,
                }
            )
    return pd.DataFrame(rows)


def build_hist(series: pd.Series, bins: list[float], metric_name: str) -> pd.DataFrame:
    binned = pd.cut(series, bins=bins, right=False, include_lowest=True)
    counts = binned.value_counts().sort_index()
    out = counts.rename_axis("bin").reset_index(name="count")
    out["fraction"] = out["count"] / max(int(series.shape[0]), 1)
    out["metric"] = metric_name
    return out[["metric", "bin", "count", "fraction"]]


def main() -> None:
    args = parse_args()
    reads = Path(args.reads)
    ref_mmi = Path(args.ref_mmi)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    if not reads.exists():
        raise FileNotFoundError(f"Reads file not found: {reads}")
    if not ref_mmi.exists():
        raise FileNotFoundError(f"Reference index not found: {ref_mmi}")

    ensure_tool("minimap2")
    ensure_tool("gzip")

    paf_path = Path(f"{out_prefix}.alignment.paf")
    paf_gz_path = Path(f"{out_prefix}.alignment.paf.gz")

    run_minimap2_paf(reads, ref_mmi, args.threads, paf_path)

    with paf_path.open("rb") as src, gzip.open(paf_gz_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    if not args.keep_paf:
        paf_path.unlink(missing_ok=True)

    read_lengths = dict(iter_fastq_lengths(reads))
    per_read_base = pd.DataFrame(
        {
            "read_id": list(read_lengths.keys()),
            "read_length": list(read_lengths.values()),
        }
    )

    cols = [
        "read_id",
        "query_len",
        "query_start",
        "query_end",
        "strand",
        "target",
        "target_len",
        "target_start",
        "target_end",
        "nmatch",
        "aln_block_len",
        "mapq",
    ]
    if paf_gz_path.stat().st_size > 0:
        paf_df = pd.read_csv(paf_gz_path, sep="\t", header=None, usecols=range(12), names=cols, compression="gzip")
    else:
        paf_df = pd.DataFrame(columns=cols)

    if not paf_df.empty:
        paf_df["nmatch"] = paf_df["nmatch"].astype(float)
        paf_df["aln_block_len"] = paf_df["aln_block_len"].astype(float)
        paf_df["query_len"] = paf_df["query_len"].astype(float)
        paf_df["query_start"] = paf_df["query_start"].astype(float)
        paf_df["query_end"] = paf_df["query_end"].astype(float)
        paf_df["aligned_query_bases"] = (paf_df["query_end"] - paf_df["query_start"]).clip(lower=0)
        paf_df["identity"] = (paf_df["nmatch"] / paf_df["aln_block_len"]).fillna(0.0)
        paf_df["query_cov"] = (paf_df["aligned_query_bases"] / paf_df["query_len"]).fillna(0.0)

        best = select_best_alignment(paf_df)

        target_meta = best["target"].map(parse_gencode_target_name).apply(pd.Series)
        best = pd.concat([best, target_meta], axis=1)
        best = best.rename(
            columns={
                "target": "best_target",
                "query_cov": "best_query_cov",
                "identity": "best_identity",
                "nmatch": "best_nmatch",
                "aln_block_len": "best_aln_block_len",
                "mapq": "best_mapq",
            }
        )
        keep_cols = [
            "read_id",
            "best_target",
            "best_query_cov",
            "best_identity",
            "best_nmatch",
            "best_aln_block_len",
            "best_mapq",
            "transcript_id",
            "gene_id",
            "transcript_name",
            "gene_name",
            "biotype",
        ]
        best = best[keep_cols]
    else:
        best = pd.DataFrame(columns=[
            "read_id",
            "best_target",
            "best_query_cov",
            "best_identity",
            "best_nmatch",
            "best_aln_block_len",
            "best_mapq",
            "transcript_id",
            "gene_id",
            "transcript_name",
            "gene_name",
            "biotype",
        ])

    per_read = per_read_base.merge(best, on="read_id", how="left")
    per_read["best_query_cov"] = per_read["best_query_cov"].fillna(0.0)
    per_read["best_identity"] = per_read["best_identity"].fillna(0.0)
    per_read["mapped_to_human"] = per_read["best_target"].notna()
    per_read["explained"] = (
        (per_read["best_query_cov"] >= args.min_query_cov)
        & (per_read["best_identity"] >= args.min_identity)
    )

    per_read_path = Path(f"{out_prefix}.per_read.tsv")
    per_read.to_csv(per_read_path, sep="\t", index=False)

    total_reads = len(per_read)
    explained_reads = int(per_read["explained"].sum())
    mapped_reads = int(per_read["mapped_to_human"].sum())
    total_bases = int(per_read["read_length"].sum())
    explained_bases = int(per_read.loc[per_read["explained"], "read_length"].sum())

    summary_rows = [
        {"metric": "total_reads", "value": total_reads},
        {"metric": "explained_reads", "value": explained_reads},
        {"metric": "reads_explained_fraction", "value": explained_reads / max(total_reads, 1)},
        {"metric": "mapped_to_human_reads", "value": mapped_reads},
        {"metric": "mapped_to_human_fraction", "value": mapped_reads / max(total_reads, 1)},
        {"metric": "total_bases", "value": total_bases},
        {"metric": "explained_bases", "value": explained_bases},
        {"metric": "bases_explained_fraction", "value": explained_bases / max(total_bases, 1)},
        {"metric": "min_query_cov", "value": args.min_query_cov},
        {"metric": "min_identity", "value": args.min_identity},
    ]

    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        summary_rows.append(
            {
                "metric": f"best_identity_quantile_{q}",
                "value": float(per_read["best_identity"].quantile(q)),
            }
        )
        summary_rows.append(
            {
                "metric": f"best_query_cov_quantile_{q}",
                "value": float(per_read["best_query_cov"].quantile(q)),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary_path = Path(f"{out_prefix}.summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    bins = [0, 100, 200, 400, 800, 1200, 2000, 5000, 1000000000]
    labels = ["0-99", "100-199", "200-399", "400-799", "800-1199", "1200-1999", "2000-4999", "5000+"]
    per_read["length_bin"] = pd.cut(per_read["read_length"], bins=bins, labels=labels, right=False)
    length_binned = (
        per_read.groupby("length_bin", observed=True)
        .agg(
            total_reads=("read_id", "count"),
            explained_reads=("explained", "sum"),
            mapped_reads=("mapped_to_human", "sum"),
            mean_best_identity=("best_identity", "mean"),
            mean_best_query_cov=("best_query_cov", "mean"),
        )
        .reset_index()
    )
    length_binned["explained_fraction"] = length_binned["explained_reads"] / length_binned["total_reads"]
    length_binned["mapped_fraction"] = length_binned["mapped_reads"] / length_binned["total_reads"]
    length_binned_path = Path(f"{out_prefix}.length_binned.tsv")
    length_binned.to_csv(length_binned_path, sep="\t", index=False)

    mapped_only = per_read[per_read["mapped_to_human"]].copy()
    if not mapped_only.empty:
        top_targets = (
            mapped_only.groupby(
                ["transcript_id", "gene_id", "gene_name", "transcript_name", "biotype"],
                dropna=False,
            )
            .agg(
                mapped_reads=("read_id", "count"),
                explained_reads=("explained", "sum"),
                mean_best_identity=("best_identity", "mean"),
                mean_best_query_cov=("best_query_cov", "mean"),
            )
            .reset_index()
            .sort_values(["explained_reads", "mapped_reads", "mean_best_query_cov"], ascending=[False, False, False])
        )
    else:
        top_targets = pd.DataFrame(
            columns=[
                "transcript_id",
                "gene_id",
                "gene_name",
                "transcript_name",
                "biotype",
                "mapped_reads",
                "explained_reads",
                "mean_best_identity",
                "mean_best_query_cov",
            ]
        )
    top_targets_path = Path(f"{out_prefix}.top_targets.tsv")
    top_targets.to_csv(top_targets_path, sep="\t", index=False)

    grid = build_threshold_grid(per_read)
    grid_path = Path(f"{out_prefix}.threshold_grid.tsv")
    grid.to_csv(grid_path, sep="\t", index=False)

    identity_hist = build_hist(per_read["best_identity"], [x / 100.0 for x in range(0, 105, 5)], "best_identity")
    identity_hist_path = Path(f"{out_prefix}.identity_hist.tsv")
    identity_hist.to_csv(identity_hist_path, sep="\t", index=False)

    qcov_hist = build_hist(per_read["best_query_cov"], [x / 100.0 for x in range(0, 105, 5)], "best_query_cov")
    qcov_hist_path = Path(f"{out_prefix}.query_cov_hist.tsv")
    qcov_hist.to_csv(qcov_hist_path, sep="\t", index=False)

    print(f"Wrote: {per_read_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {length_binned_path}")
    print(f"Wrote: {top_targets_path}")
    print(f"Wrote: {grid_path}")
    print(f"Wrote: {identity_hist_path}")
    print(f"Wrote: {qcov_hist_path}")
    print(f"Wrote: {paf_gz_path}")


if __name__ == "__main__":
    main()

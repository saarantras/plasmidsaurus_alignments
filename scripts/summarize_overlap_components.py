#!/usr/bin/env python3
"""Summarize overlap-graph component structure from PAF alignments."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--paf", required=True, help="PAF file from all-vs-all overlap mapping")
    p.add_argument("--read-ids", default="", help="Optional newline-delimited read IDs to force full universe")
    p.add_argument("--min-overlap", type=int, default=150, help="Minimum alignment block length")
    p.add_argument("--min-identity", type=float, default=0.85, help="Minimum identity (nmatch/aln_block_len)")
    p.add_argument("--out-prefix", required=True, help="Output prefix")
    return p.parse_args()


def find(parent: dict[str, str], x: str) -> str:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(parent: dict[str, str], size: dict[str, int], a: str, b: str) -> None:
    ra = find(parent, a)
    rb = find(parent, b)
    if ra == rb:
        return
    if size[ra] < size[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    size[ra] += size[rb]


def main() -> None:
    args = parse_args()
    paf_path = Path(args.paf)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    if not paf_path.exists():
        raise FileNotFoundError(f"PAF not found: {paf_path}")

    all_ids: set[str] = set()
    if args.read_ids:
        rid_path = Path(args.read_ids)
        if not rid_path.exists():
            raise FileNotFoundError(f"Read ID file not found: {rid_path}")
        all_ids = {line.strip() for line in rid_path.open() if line.strip()}

    parent: dict[str, str] = {}
    size: dict[str, int] = {}
    seen: set[str] = set()
    edges_kept = 0
    edges_total = 0

    # Initialize read universe if explicitly provided.
    for rid in all_ids:
        parent[rid] = rid
        size[rid] = 1

    with paf_path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q = f[0]
            t = f[5]

            if all_ids and (q not in all_ids or t not in all_ids):
                continue

            seen.add(q)
            seen.add(t)

            if q not in parent:
                parent[q] = q
                size[q] = 1
            if t not in parent:
                parent[t] = t
                size[t] = 1

            if q == t:
                continue

            edges_total += 1
            nmatch = int(f[9])
            block_len = int(f[10])
            if block_len == 0:
                continue
            ident = nmatch / block_len
            if block_len >= args.min_overlap and ident >= args.min_identity:
                union(parent, size, q, t)
                edges_kept += 1

    universe = all_ids if all_ids else set(parent.keys())
    comp = Counter(find(parent, rid) for rid in universe)
    sizes = sorted(comp.values(), reverse=True)
    if not sizes:
        sizes = [0]

    summary = pd.DataFrame(
        [
            {"metric": "reads_total", "value": len(universe)},
            {"metric": "reads_seen_in_paf", "value": len(seen)},
            {"metric": "reads_no_any_overlap_record", "value": len(universe) - len(seen)},
            {"metric": "min_overlap", "value": args.min_overlap},
            {"metric": "min_identity", "value": args.min_identity},
            {"metric": "edges_total_nonself", "value": edges_total},
            {"metric": "edges_kept_filtered", "value": edges_kept},
            {"metric": "n_components", "value": len(sizes)},
            {"metric": "largest_component_reads", "value": sizes[0]},
            {"metric": "largest_component_fraction", "value": (sizes[0] / len(universe)) if universe else 0.0},
            {"metric": "singletons", "value": sum(1 for s in sizes if s == 1)},
            {"metric": "components_ge_10", "value": sum(1 for s in sizes if s >= 10)},
        ]
    )
    summary_path = Path(f"{out_prefix}.summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    comp_df = pd.DataFrame({"component_size": sizes})
    comp_df["component_rank"] = range(1, len(comp_df) + 1)
    comp_df = comp_df[["component_rank", "component_size"]]
    comp_path = Path(f"{out_prefix}.component_sizes.tsv")
    comp_df.to_csv(comp_path, sep="\t", index=False)

    print(f"Wrote: {summary_path}")
    print(f"Wrote: {comp_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Convert a GenBank file into FASTA and BED feature tracks."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO


ID_KEYS = ("label", "gene", "locus_tag", "note", "product")


def pick_feature_id(feature, fallback: str) -> str:
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    for key in ID_KEYS:
        values = qualifiers.get(key)
        if values:
            value = values[0] if isinstance(values, list) else values
            cleaned = str(value).replace("\t", " ").replace("\n", " ").strip()
            if cleaned:
                return cleaned
    return fallback


def strand_symbol(location) -> str:
    strand = getattr(location, "strand", None)
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genbank", required=True, help="Input GenBank file")
    parser.add_argument("--fasta", required=True, help="Output FASTA path")
    parser.add_argument("--bed", required=True, help="Output BED path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genbank_path = Path(args.genbank)
    fasta_path = Path(args.fasta)
    bed_path = Path(args.bed)

    record = SeqIO.read(genbank_path, "genbank")

    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    bed_path.parent.mkdir(parents=True, exist_ok=True)

    SeqIO.write(record, fasta_path, "fasta")

    rows_written = 0
    with bed_path.open("w", encoding="utf-8") as bed_out:
        for feat_idx, feature in enumerate(record.features, start=1):
            if feature.type == "source":
                continue

            feature_id = pick_feature_id(feature, f"{feature.type}_{feat_idx}")
            parts = getattr(feature.location, "parts", [feature.location])
            strand = strand_symbol(feature.location)

            for part_idx, part in enumerate(parts, start=1):
                start = int(part.start)
                end = int(part.end)
                if end <= start:
                    continue

                part_feature_id = feature_id
                if len(parts) > 1:
                    part_feature_id = f"{feature_id}_part{part_idx}"

                fields = [
                    record.id,
                    str(start),
                    str(end),
                    feature.type,
                    part_feature_id,
                    strand,
                ]
                bed_out.write("\t".join(fields) + "\n")
                rows_written += 1

    print(f"Wrote FASTA: {fasta_path}")
    print(f"Wrote BED features: {bed_path} ({rows_written} intervals)")


if __name__ == "__main__":
    main()

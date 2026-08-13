#!/usr/bin/env python3
"""Export per-read observations for PASS SNVs from both orientations.

Rows from the two orientations are deliberately kept separate; do not sum them.
"""

import argparse
import csv
from pathlib import Path

import pysam


MTDNA_LENGTH = 16_299
SPHI_FIRST_STANDARD_POSITION = 10_759
NCOI_FIRST_STANDARD_POSITION = 9_216


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pass-calls", type=Path, required=True)
    p.add_argument("--sphi-bam", type=Path, required=True)
    p.add_argument("--ncoi-bam", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--sphi-contig", default="NC_005089.1_SphI_rotated")
    p.add_argument("--ncoi-contig", default="NC_005089.1_NcoI_rotated")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=20)
    return p.parse_args()


def standard_to_rotated(position, first_standard_position):
    return ((position - first_standard_position) % MTDNA_LENGTH) + 1


def load_calls(path):
    with path.open(newline="") as handle:
        calls = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_id", "standard_position", "ref", "alt", "selected_orientation"}
    if calls and not required.issubset(calls[0]):
        raise ValueError("PASS table is missing required columns")
    return calls


def observations(bam_path, contig, orientation, first_standard_position, call, min_mapq, min_baseq):
    standard_position = int(call["standard_position"])
    rotated_position = standard_to_rotated(standard_position, first_standard_position)
    rows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for column in bam.pileup(contig, rotated_position - 1, rotated_position,
                                 truncate=True, stepper="samtools",
                                 min_mapping_quality=min_mapq, min_base_quality=0,
                                 ignore_overlaps=False, ignore_orphans=False, max_depth=100_000):
            if column.reference_pos != rotated_position - 1:
                continue
            for pileup_read in column.pileups:
                read = pileup_read.alignment
                if (read.is_unmapped or read.is_secondary or read.is_supplementary
                        or read.is_duplicate or read.is_qcfail or read.mapping_quality < min_mapq
                        or pileup_read.is_del or pileup_read.is_refskip
                        or pileup_read.query_position is None or read.query_sequence is None
                        or read.query_qualities is None):
                    continue
                qpos = pileup_read.query_position
                baseq = read.query_qualities[qpos]
                if baseq < min_baseq:
                    continue
                base = read.query_sequence[qpos].upper()
                rows.append({
                    "sample_id": call["sample_id"], "standard_position": standard_position,
                    "ref": call["ref"], "alt": call["alt"], "orientation": orientation,
                    "is_selected_orientation": str(orientation == call["selected_orientation"]),
                    "read_id": read.query_name, "observed_base": base,
                    "supports_ref": str(base == call["ref"]),
                    "supports_alt": str(base == call["alt"]),
                    "strand": "reverse" if read.is_reverse else "forward",
                    "mapq": read.mapping_quality, "baseq": baseq,
                    "query_position_0based": qpos,
                    "query_length": read.query_length or "",
                    "cigar": read.cigarstring or "",
                    "alignment_start_1based": read.reference_start + 1,
                    "alignment_end_1based": read.reference_end,
                })
    return rows


def main():
    args = parse_args()
    calls = load_calls(args.pass_calls)
    fields = ["sample_id", "standard_position", "ref", "alt", "orientation",
              "is_selected_orientation", "read_id", "observed_base", "supports_ref",
              "supports_alt", "strand", "mapq", "baseq", "query_position_0based",
              "query_length", "cigar", "alignment_start_1based", "alignment_end_1based"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for call in calls:
            writer.writerows(observations(args.sphi_bam, args.sphi_contig, "SphI",
                                          SPHI_FIRST_STANDARD_POSITION, call,
                                          args.min_mapq, args.min_baseq))
            writer.writerows(observations(args.ncoi_bam, args.ncoi_contig, "NcoI",
                                          NCOI_FIRST_STANDARD_POSITION, call,
                                          args.min_mapq, args.min_baseq))
    print(f"Exported read-level support for {len(calls)} PASS SNV(s) to {args.output}")


if __name__ == "__main__":
    main()

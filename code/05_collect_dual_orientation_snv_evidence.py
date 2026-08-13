#!/usr/bin/env python3
"""Collect CIGAR-aware SNV evidence from SphI- and NcoI-oriented mtDNA BAMs.

Each molecule is observed in both orientations, but counts are kept separate.
The selected evidence is SphI by default and NcoI near the SphI boundary.
"""

import argparse
import csv
from pathlib import Path

import pysam


MTDNA_LENGTH = 16_299
SPHI_FIRST_STANDARD_POSITION = 10_759
NCOI_FIRST_STANDARD_POSITION = 9_216
SPHI_BOUNDARIES = (10_758, 10_759)
NCOI_BOUNDARIES = (9_215, 9_216)
BASES = ("A", "C", "G", "T")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--sphi-bam", type=Path, required=True)
    p.add_argument("--ncoi-bam", type=Path, required=True)
    p.add_argument("--sphi-reference", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--sphi-contig", default="NC_005089.1_SphI_rotated")
    p.add_argument("--ncoi-contig", default="NC_005089.1_NcoI_rotated")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=20)
    p.add_argument("--boundary-window", type=int, default=100)
    p.add_argument("--max-depth", type=int, default=100_000)
    return p.parse_args()


def rotate_to_standard(position, first_standard_position):
    return ((position + first_standard_position - 2) % MTDNA_LENGTH) + 1


def standard_to_rotated(position, first_standard_position):
    return ((position - first_standard_position) % MTDNA_LENGTH) + 1


def circular_distance(position, boundary):
    difference = abs(position - boundary)
    return min(difference, MTDNA_LENGTH - difference)


def near_boundary(position, boundaries, window):
    return any(circular_distance(position, x) <= window for x in boundaries)


def empty_record():
    record = {"depth": 0}
    for base in BASES:
        record[base] = 0
        record[f"{base}_forward"] = 0
        record[f"{base}_reverse"] = 0
        record[f"{base}_baseq_sum"] = 0
    return record


def validate_bam(bam_path, contig):
    if not bam_path.exists():
        raise FileNotFoundError(bam_path)
    if not Path(str(bam_path) + ".bai").exists() and not bam_path.with_suffix(".bai").exists():
        raise FileNotFoundError(f"BAM index not found for {bam_path}")
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if contig not in bam.references:
            raise ValueError(f"{contig!r} is absent from {bam_path}")
        if bam.get_reference_length(contig) != MTDNA_LENGTH:
            raise ValueError(f"{contig!r} is not {MTDNA_LENGTH:,} bp in {bam_path}")


def count_orientation(bam_path, contig, first_standard_position, min_mapq, min_baseq, max_depth):
    counts = {position: empty_record() for position in range(1, MTDNA_LENGTH + 1)}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for column in bam.pileup(
            contig, 0, MTDNA_LENGTH, truncate=True, stepper="samtools",
            min_mapping_quality=min_mapq, min_base_quality=0,
            ignore_overlaps=False, ignore_orphans=False, max_depth=max_depth,
        ):
            standard_position = rotate_to_standard(column.reference_pos + 1, first_standard_position)
            record = counts[standard_position]
            for pileup_read in column.pileups:
                read = pileup_read.alignment
                if (
                    read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.is_qcfail or read.mapping_quality < min_mapq
                    or pileup_read.is_del or pileup_read.is_refskip
                    or pileup_read.query_position is None or read.query_sequence is None
                    or read.query_qualities is None
                ):
                    continue
                qpos = pileup_read.query_position
                baseq = read.query_qualities[qpos]
                if baseq < min_baseq:
                    continue
                base = read.query_sequence[qpos].upper()
                if base not in BASES:
                    continue
                strand = "reverse" if read.is_reverse else "forward"
                record["depth"] += 1
                record[base] += 1
                record[f"{base}_{strand}"] += 1
                record[f"{base}_baseq_sum"] += baseq
    return counts


def prefixed(record, prefix):
    output = {f"{prefix}_depth": record["depth"]}
    for base in BASES:
        output[f"{prefix}_{base}"] = record[base]
        output[f"{prefix}_{base}_forward"] = record[f"{base}_forward"]
        output[f"{prefix}_{base}_reverse"] = record[f"{base}_reverse"]
        output[f"{prefix}_{base}_mean_baseq"] = (
            record[f"{base}_baseq_sum"] / record[base] if record[base] else 0.0
        )
    return output


def read_standard_reference(path):
    if not path.exists():
        raise FileNotFoundError(path)
    with pysam.FastaFile(str(path)) as fasta:
        if len(fasta.references) != 1:
            raise ValueError("--sphi-reference must contain exactly one mtDNA contig")
        sequence = fasta.fetch(fasta.references[0]).upper()
    if len(sequence) != MTDNA_LENGTH:
        raise ValueError(f"Expected {MTDNA_LENGTH:,} bp; found {len(sequence):,}")
    standard = [None] * MTDNA_LENGTH
    for rotated_position, base in enumerate(sequence, start=1):
        standard_position = rotate_to_standard(rotated_position, SPHI_FIRST_STANDARD_POSITION)
        standard[standard_position - 1] = base
    return "".join(standard)


def main():
    args = parse_args()
    for value, name in ((args.min_mapq, "min-mapq"), (args.min_baseq, "min-baseq"),
                        (args.boundary_window, "boundary-window"), (args.max_depth, "max-depth")):
        if value < 0:
            raise ValueError(f"--{name} must be nonnegative")
    validate_bam(args.sphi_bam, args.sphi_contig)
    validate_bam(args.ncoi_bam, args.ncoi_contig)
    reference = read_standard_reference(args.sphi_reference)

    print(f"{args.sample}: counting SphI orientation", flush=True)
    sphi = count_orientation(args.sphi_bam, args.sphi_contig, SPHI_FIRST_STANDARD_POSITION,
                             args.min_mapq, args.min_baseq, args.max_depth)
    print(f"{args.sample}: counting NcoI orientation", flush=True)
    ncoi = count_orientation(args.ncoi_bam, args.ncoi_contig, NCOI_FIRST_STANDARD_POSITION,
                             args.min_mapq, args.min_baseq, args.max_depth)

    metadata = ["sample_id", "standard_position", "ref", "context_11bp", "selected_orientation",
                "validation_orientation", "near_sphi_boundary", "near_ncoi_boundary",
                "min_mapq", "min_baseq"]
    orientation_fields = []
    for prefix in ("sphi", "ncoi"):
        orientation_fields.append(f"{prefix}_depth")
        for base in BASES:
            orientation_fields.extend([f"{prefix}_{base}", f"{prefix}_{base}_forward",
                                       f"{prefix}_{base}_reverse", f"{prefix}_{base}_mean_baseq"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata + orientation_fields, delimiter="\t")
        writer.writeheader()
        for position in range(1, MTDNA_LENGTH + 1):
            near_sphi = near_boundary(position, SPHI_BOUNDARIES, args.boundary_window)
            near_ncoi = near_boundary(position, NCOI_BOUNDARIES, args.boundary_window)
            selected = "NcoI" if near_sphi else "SphI"
            row = {
                "sample_id": args.sample, "standard_position": position,
                "ref": reference[position - 1], "selected_orientation": selected,
                "context_11bp": "".join(reference[(position - 1 + offset) % MTDNA_LENGTH]
                                          for offset in range(-5, 6)),
                "validation_orientation": "SphI" if selected == "NcoI" else "NcoI",
                "near_sphi_boundary": str(near_sphi), "near_ncoi_boundary": str(near_ncoi),
                "min_mapq": args.min_mapq, "min_baseq": args.min_baseq,
            }
            row.update(prefixed(sphi[position], "sphi"))
            row.update(prefixed(ncoi[position], "ncoi"))
            writer.writerow(row)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

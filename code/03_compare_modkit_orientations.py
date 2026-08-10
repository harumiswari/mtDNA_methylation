#!/usr/bin/env python3

import csv
import gzip
from collections import defaultdict
from pathlib import Path

SAMPLES = [f"XX_{number:02d}" for number in range(1, 9)]

INPUT_DIR = Path("/home/harumiswari/Documents/mtdna/XX_methylation_modkit_results")
OUTPUT_DIR = Path("/home/harumiswari/Documents/mtdna/XX_methylation_comparison_results")

MTDNA_LENGTH = 16_299
SPHI_FIRST_STANDARD_POSITION = 10_759
NCOI_FIRST_SPHI_POSITION = 14_757
NCOI_BOUNDARY_STANDARD_POSITION = 9_216
BOUNDARY_WINDOW = 100


def rotated_to_standard(position, first_standard_position):
    return ((position + first_standard_position - 2) % MTDNA_LENGTH) + 1


def standard_to_sphi(position):
    return ((position - SPHI_FIRST_STANDARD_POSITION) % MTDNA_LENGTH) + 1


def standard_to_ncoi(position):
    sphi_position = standard_to_sphi(position)
    return ((sphi_position - NCOI_FIRST_SPHI_POSITION) % MTDNA_LENGTH) + 1


def ncoi_to_standard(position):
    sphi_position = ((position + NCOI_FIRST_SPHI_POSITION - 2) % MTDNA_LENGTH) + 1
    return rotated_to_standard(sphi_position, SPHI_FIRST_STANDARD_POSITION)


def read_modkit_file(path, orientation):
    records = defaultdict(dict)

    with gzip.open(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 18:
                raise ValueError(f"Unexpected modkit row with {len(row)} columns: {path}")

            rotated_position = int(row[1]) + 1
            if orientation == "SphI":
                standard_position = rotated_to_standard(
                    rotated_position,
                    SPHI_FIRST_STANDARD_POSITION,
                )
            elif orientation == "NcoI":
                standard_position = ncoi_to_standard(rotated_position)
            else:
                raise ValueError(f"Unknown orientation: {orientation}")

            modification = row[3]
            if modification == "m":
                label = "5mC"
            elif modification == "h":
                label = "5hmC"
            else:
                continue

            coverage = int(row[9])
            fraction_modified = float(row[10])
            modified_count = int(row[11])
            canonical_count = int(row[12])
            other_count = int(row[13])
            delete_count = int(row[14])
            fail_count = int(row[15])
            diff_count = int(row[16])

            records[standard_position][label] = {
                "rotated_position": rotated_position,
                "coverage": coverage,
                "fraction_modified": fraction_modified,
                "modified_count": modified_count,
                "canonical_count": canonical_count,
                "other_count": other_count,
                "delete_count": delete_count,
                "fail_count": fail_count,
                "diff_count": diff_count,
            }

    return records


def write_sample_table(sample_id, sphi, ncoi):
    positions = sorted(set(sphi) | set(ncoi))
    output = OUTPUT_DIR / f"{sample_id}.SphI_NcoI.comparison.tsv"

    fields = [
        "sample_id",
        "standard_position",
        "sphI_position",
        "ncoI_position",
        "near_sphI_boundary",
        "near_ncoI_boundary",
        "sphI_5mC_coverage",
        "sphI_5mC_modified_count",
        "sphI_5mC_fraction",
        "sphI_5hmC_modified_count",
        "sphI_5hmC_fraction",
        "sphI_total_modified_fraction",
        "ncoI_5mC_coverage",
        "ncoI_5mC_modified_count",
        "ncoI_5mC_fraction",
        "ncoI_5hmC_modified_count",
        "ncoI_5hmC_fraction",
        "ncoI_total_modified_fraction",
        "status",
    ]

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for position in positions:
            sph_data = sphi.get(position, {})
            nco_data = ncoi.get(position, {})
            sph_m = sph_data.get("5mC", {})
            sph_h = sph_data.get("5hmC", {})
            nco_m = nco_data.get("5mC", {})
            nco_h = nco_data.get("5hmC", {})

            sph_coverage = sph_m.get("coverage", sph_h.get("coverage", 0))
            nco_coverage = nco_m.get("coverage", nco_h.get("coverage", 0))
            sph_total = (
                (sph_m.get("modified_count", 0) + sph_h.get("modified_count", 0))
                / sph_coverage
                if sph_coverage else 0.0
            )
            nco_total = (
                (nco_m.get("modified_count", 0) + nco_h.get("modified_count", 0))
                / nco_coverage
                if nco_coverage else 0.0
            )

            near_sphI = abs(position - 10_759) <= BOUNDARY_WINDOW or abs(position - 16_299) <= BOUNDARY_WINDOW
            near_ncoI = abs(position - NCOI_BOUNDARY_STANDARD_POSITION) <= BOUNDARY_WINDOW

            if near_sphI or near_ncoI:
                status = "boundary_review"
            elif sph_coverage == 0 or nco_coverage == 0:
                status = "one_orientation_missing"
            else:
                status = "compare"

            writer.writerow(
                {
                    "sample_id": sample_id,
                    "standard_position": position,
                    "sphI_position": ((position - SPHI_FIRST_STANDARD_POSITION) % MTDNA_LENGTH) + 1,
                    "ncoI_position": standard_to_ncoi(position),
                    "near_sphI_boundary": str(near_sphI),
                    "near_ncoI_boundary": str(near_ncoI),
                    "sphI_5mC_coverage": sph_coverage,
                    "sphI_5mC_modified_count": sph_m.get("modified_count", 0),
                    "sphI_5mC_fraction": sph_m.get("fraction_modified", 0.0),
                    "sphI_5hmC_modified_count": sph_h.get("modified_count", 0),
                    "sphI_5hmC_fraction": sph_h.get("fraction_modified", 0.0),
                    "sphI_total_modified_fraction": sph_total,
                    "ncoI_5mC_coverage": nco_coverage,
                    "ncoI_5mC_modified_count": nco_m.get("modified_count", 0),
                    "ncoI_5mC_fraction": nco_m.get("fraction_modified", 0.0),
                    "ncoI_5hmC_modified_count": nco_h.get("modified_count", 0),
                    "ncoI_5hmC_fraction": nco_h.get("fraction_modified", 0.0),
                    "ncoI_total_modified_fraction": nco_total,
                    "status": status,
                }
            )

    return output


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for sample_id in SAMPLES:
        sphi_path = INPUT_DIR / f"{sample_id}.SphI.5mC_5hmC.CpG.bed.gz"
        ncoi_path = INPUT_DIR / f"{sample_id}.NcoI.5mC_5hmC.CpG.bed.gz"

        if not sphi_path.exists() or not ncoi_path.exists():
            print(f"Skipping {sample_id}: modkit output is missing")
            continue

        sphi = read_modkit_file(sphi_path, "SphI")
        ncoi = read_modkit_file(ncoi_path, "NcoI")
        output = write_sample_table(sample_id, sphi, ncoi)

        summary_rows.append(
            {
                "sample_id": sample_id,
                "positions_sphI": len(sphi),
                "positions_ncoI": len(ncoi),
                "positions_both": len(set(sphi) & set(ncoi)),
                "output": output,
            }
        )
        print(f"Wrote {output}")

    summary_path = OUTPUT_DIR / "XX_replicate_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "positions_sphI", "positions_ncoI", "positions_both", "output"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()

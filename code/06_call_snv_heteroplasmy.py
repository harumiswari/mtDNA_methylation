#!/usr/bin/env python3
"""Call and QC mtDNA SNV heteroplasmy from dual-orientation evidence."""

import argparse
import csv
import json
import math
from pathlib import Path


BASES = ("A", "C", "G", "T")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--output-prefix", type=Path, required=True)
    p.add_argument("--min-depth", type=int, default=100)
    p.add_argument("--min-alt-count", type=int, default=5)
    p.add_argument("--min-vaf", type=float, default=0.01, help="Fraction; 0.01 means 1%%")
    p.add_argument("--max-vaf", type=float, default=0.95, help="Exclude near-homoplasmic/reference differences")
    p.add_argument("--min-alt-per-strand", type=int, default=2)
    p.add_argument("--min-alt-mean-baseq", type=float, default=20.0)
    p.add_argument("--strand-bias-p", type=float, default=0.001)
    p.add_argument("--min-validation-depth", type=int, default=100)
    p.add_argument("--max-orientation-vaf-diff", type=float, default=0.02)
    p.add_argument("--max-orientation-relative-diff", type=float, default=0.25)
    p.add_argument("--max-homopolymer-run", type=int, default=6,
                   help="Fail reference homopolymer runs longer than this; 0 disables")
    p.add_argument("--allow-missing-validation", action="store_true")
    return p.parse_args()


def require_columns(fieldnames):
    required = {"sample_id", "standard_position", "ref", "context_11bp",
                "selected_orientation", "near_sphi_boundary", "near_ncoi_boundary"}
    for prefix in ("sphi", "ncoi"):
        required.add(f"{prefix}_depth")
        for base in BASES:
            required.update({f"{prefix}_{base}", f"{prefix}_{base}_forward",
                             f"{prefix}_{base}_reverse", f"{prefix}_{base}_mean_baseq"})
    missing = sorted(required - set(fieldnames or []))
    if missing:
        raise ValueError("Evidence table is missing: " + ", ".join(missing))


def as_bool(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_two_sided(a, b, c, d):
    """Two-sided Fisher exact p for [[a,b],[c,d]], without scipy."""
    row1, row2, col1 = a + b, c + d, a + c
    total = row1 + row2
    if total == 0:
        return 1.0
    low, high = max(0, col1 - row2), min(row1, col1)

    def probability(x):
        return math.exp(log_choose(col1, x) + log_choose(total - col1, row1 - x)
                        - log_choose(total, row1))

    observed = probability(a)
    return min(1.0, sum(probability(x) for x in range(low, high + 1)
                        if probability(x) <= observed + 1e-12))


def homopolymer_run(context):
    if not context:
        return 0
    center = len(context) // 2
    base = context[center]
    left = center
    right = center
    while left > 0 and context[left - 1] == base:
        left -= 1
    while right + 1 < len(context) and context[right + 1] == base:
        right += 1
    return right - left + 1


def orientation_values(row, orientation, alt):
    prefix = orientation.lower()
    depth = int(row[f"{prefix}_depth"])
    alt_count = int(row[f"{prefix}_{alt}"])
    return {
        "depth": depth,
        "alt_count": alt_count,
        "vaf": alt_count / depth if depth else 0.0,
        "alt_forward": int(row[f"{prefix}_{alt}_forward"]),
        "alt_reverse": int(row[f"{prefix}_{alt}_reverse"]),
        "alt_mean_baseq": float(row[f"{prefix}_{alt}_mean_baseq"]),
    }


def call_row(row, args):
    ref = row["ref"].upper()
    if ref not in BASES:
        raise ValueError(f"Non-ACGT reference base at position {row['standard_position']}: {ref}")
    selected = row["selected_orientation"]
    primary_prefix = selected.lower()
    validation = "NcoI" if selected == "SphI" else "SphI"
    depth = int(row[f"{primary_prefix}_depth"])
    alt = max((base for base in BASES if base != ref),
              key=lambda base: int(row[f"{primary_prefix}_{base}"]))
    primary = orientation_values(row, selected, alt)
    secondary = orientation_values(row, validation, alt)
    ref_forward = int(row[f"{primary_prefix}_{ref}_forward"])
    ref_reverse = int(row[f"{primary_prefix}_{ref}_reverse"])
    strand_p = fisher_exact_two_sided(primary["alt_forward"], primary["alt_reverse"],
                                      ref_forward, ref_reverse)
    near_boundary = as_bool(row["near_sphi_boundary"]) or as_bool(row["near_ncoi_boundary"])
    vaf_difference = abs(primary["vaf"] - secondary["vaf"])
    allowed_difference = max(args.max_orientation_vaf_diff,
                             args.max_orientation_relative_diff * primary["vaf"])
    run = homopolymer_run(row["context_11bp"].upper())

    filters = []
    if primary["depth"] < args.min_depth:
        filters.append("LOW_DEPTH")
    if primary["alt_count"] < args.min_alt_count:
        filters.append("LOW_ALT_COUNT")
    if primary["vaf"] < args.min_vaf:
        filters.append("LOW_VAF")
    if primary["vaf"] > args.max_vaf:
        filters.append("NEAR_HOMOPLASMIC")
    if primary["alt_forward"] < args.min_alt_per_strand:
        filters.append("LOW_ALT_FORWARD")
    if primary["alt_reverse"] < args.min_alt_per_strand:
        filters.append("LOW_ALT_REVERSE")
    if primary["alt_mean_baseq"] < args.min_alt_mean_baseq:
        filters.append("LOW_ALT_MEAN_BASEQ")
    if strand_p < args.strand_bias_p:
        filters.append("STRAND_BIAS")
    if args.max_homopolymer_run and run > args.max_homopolymer_run:
        filters.append("HOMOPOLYMER_CONTEXT")
    # The alternate orientation is distorted at one of the artificial boundaries,
    # so concordance is assessed only away from both boundary windows.
    if not near_boundary:
        if secondary["depth"] < args.min_validation_depth:
            if not args.allow_missing_validation:
                filters.append("LOW_VALIDATION_DEPTH")
        elif vaf_difference > allowed_difference:
            filters.append("ORIENTATION_DISCORDANT")

    output = {
        "sample_id": row["sample_id"],
        "standard_position": int(row["standard_position"]),
        "ref": ref, "alt": alt, "change": f"{ref}>{alt}",
        "context_11bp": row["context_11bp"], "homopolymer_run": run,
        "selected_orientation": selected, "validation_orientation": validation,
        "near_sphi_boundary": row["near_sphi_boundary"],
        "near_ncoi_boundary": row["near_ncoi_boundary"],
        "depth": primary["depth"], "ref_count": int(row[f"{primary_prefix}_{ref}"]),
        "alt_count": primary["alt_count"], "vaf": primary["vaf"],
        "vaf_percent": 100.0 * primary["vaf"],
        "alt_forward": primary["alt_forward"], "alt_reverse": primary["alt_reverse"],
        "ref_forward": ref_forward, "ref_reverse": ref_reverse,
        "alt_mean_baseq": primary["alt_mean_baseq"], "strand_bias_p": strand_p,
        "validation_depth": secondary["depth"],
        "validation_alt_count": secondary["alt_count"],
        "validation_vaf": secondary["vaf"],
        "validation_vaf_percent": 100.0 * secondary["vaf"],
        "orientation_vaf_difference": vaf_difference,
        "orientation_vaf_difference_percent": 100.0 * vaf_difference,
        "filter": "PASS" if not filters else ";".join(filters),
        "pass": str(not filters),
    }
    return output


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_pass_vcf(path, rows, sample_id):
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##contig=<ID=NC_005089.1,length=16299>\n")
        handle.write('##INFO=<ID=SOURCE,Number=1,Type=String,Description="Selected alignment orientation">\n')
        handle.write('##INFO=<ID=VALAF,Number=1,Type=Float,Description="Validation-orientation alternate allele fraction">\n')
        handle.write('##INFO=<ID=OVDIFF,Number=1,Type=Float,Description="Absolute orientation AF difference">\n')
        handle.write('##INFO=<ID=SBP,Number=1,Type=Float,Description="Two-sided Fisher strand-bias p-value">\n')
        handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Selected-orientation usable depth">\n')
        handle.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reference and alternate allele depths">\n')
        handle.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele fraction">\n')
        handle.write('##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Reference and alternate forward-strand depths">\n')
        handle.write('##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Reference and alternate reverse-strand depths">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")
        for row in rows:
            info = (f"SOURCE={row['selected_orientation']};VALAF={row['validation_vaf']:.6g};"
                    f"OVDIFF={row['orientation_vaf_difference']:.6g};SBP={row['strand_bias_p']:.6g}")
            sample = (f"{row['depth']}:{row['ref_count']},{row['alt_count']}:{row['vaf']:.6g}:"
                      f"{row['ref_forward']},{row['alt_forward']}:"
                      f"{row['ref_reverse']},{row['alt_reverse']}")
            handle.write(f"NC_005089.1\t{row['standard_position']}\t.\t{row['ref']}\t{row['alt']}\t.\tPASS\t"
                         f"{info}\tDP:AD:AF:ADF:ADR\t{sample}\n")


def main():
    args = parse_args()
    if not 0 <= args.min_vaf <= args.max_vaf <= 1:
        raise ValueError("Require 0 <= --min-vaf <= --max-vaf <= 1")
    with args.evidence.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames)
        calls = [call_row(row, args) for row in reader]
    if not calls:
        raise ValueError("Evidence table contains no rows")
    fields = list(calls[0])
    screened = [x for x in calls if x["alt_count"] > 0]
    passed = [x for x in calls if x["filter"] == "PASS"]
    prefix = args.output_prefix
    all_path = Path(f"{prefix}.all_sites.tsv")
    candidate_path = Path(f"{prefix}.candidates.tsv")
    pass_path = Path(f"{prefix}.pass.tsv")
    write_tsv(all_path, calls, fields)
    write_tsv(candidate_path, screened, fields)
    write_tsv(pass_path, passed, fields)
    write_pass_vcf(Path(f"{prefix}.pass.vcf"), passed, calls[0]["sample_id"])

    filter_counts = {}
    for call in calls:
        for reason in call["filter"].split(";"):
            filter_counts[reason] = filter_counts.get(reason, 0) + 1
    summary = {
        "sample_id": calls[0]["sample_id"], "sites_evaluated": len(calls),
        "sites_with_any_alt": len(screened), "pass_snv_count": len(passed),
        "median_depth": sorted(x["depth"] for x in calls)[len(calls) // 2],
        "parameters": vars(args) | {"evidence": str(args.evidence),
                                     "output_prefix": str(args.output_prefix)},
        "filter_counts": filter_counts,
    }
    summary_path = Path(f"{prefix}.summary.json")
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2, default=str)
        handle.write("\n")
    print(f"{summary['sample_id']}: {len(passed)} PASS SNV(s); wrote {pass_path}")


if __name__ == "__main__":
    main()

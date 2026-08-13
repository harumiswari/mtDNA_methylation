#!/usr/bin/env python3
"""Create per-sample and cohort QC plots for the SNV heteroplasmy pipeline."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


MTDNA_LENGTH = 16_299
SPHI_BOUNDARY = 10_759
NCOI_BOUNDARY = 9_216


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(rows, field):
    return np.array([float(row[field]) for row in rows], dtype=float)


def mark_boundaries(axis):
    axis.axvline(SPHI_BOUNDARY, color="#d95f02", ls="--", lw=0.9, label="SphI boundary")
    axis.axvline(NCOI_BOUNDARY, color="#1b9e77", ls=":", lw=0.9, label="NcoI boundary")


def plot_sample(sample, all_rows, pass_rows, output_dir):
    positions = numeric(all_rows, "standard_position")
    depth = numeric(all_rows, "depth")
    validation_depth = numeric(all_rows, "validation_depth")
    any_alt = [row for row in all_rows if int(row["alt_count"]) > 0]

    fig, axes = plt.subplots(4, 1, figsize=(13, 11), sharex=True, constrained_layout=True)
    axes[0].plot(positions, depth, lw=0.55, color="#355f8d", label="Selected orientation")
    axes[0].plot(positions, validation_depth, lw=0.45, alpha=0.55, color="#999999",
                 label="Validation orientation")
    axes[0].set_ylabel("Usable depth")
    axes[0].legend(frameon=False, ncol=2)

    if any_alt:
        axes[1].scatter(numeric(any_alt, "standard_position"), numeric(any_alt, "vaf_percent"),
                        s=7, alpha=0.4, color="#777777", label="Any alternate call")
    if pass_rows:
        axes[1].scatter(numeric(pass_rows, "standard_position"), numeric(pass_rows, "vaf_percent"),
                        s=30, color="#b2182b", edgecolor="white", linewidth=0.4, label="PASS SNV")
    axes[1].set_ylabel("VAF (%)")
    axes[1].legend(frameon=False)

    candidate_rows = [row for row in all_rows if int(row["alt_count"]) >= 2]
    if candidate_rows:
        axes[2].scatter(numeric(candidate_rows, "standard_position"),
                        numeric(candidate_rows, "orientation_vaf_difference_percent"),
                        s=8, alpha=0.5, color="#6a3d9a")
    axes[2].set_ylabel("Orientation |ΔVAF|\n(percentage points)")

    if pass_rows:
        axes[3].scatter(numeric(pass_rows, "vaf_percent"), numeric(pass_rows, "alt_mean_baseq"),
                        s=np.maximum(20, np.sqrt(numeric(pass_rows, "alt_count")) * 8),
                        color="#e08214", alpha=0.8)
        for row in pass_rows:
            axes[3].annotate(str(row["standard_position"]),
                             (float(row["vaf_percent"]), float(row["alt_mean_baseq"])),
                             xytext=(3, 3), textcoords="offset points", fontsize=7)
    axes[3].set_xlabel("VAF (%)")
    axes[3].set_ylabel("Alternate mean base quality")
    axes[3].grid(alpha=0.2)
    axes[3].set_xlim(left=0)
    # The final panel is not genomic; do not add boundary lines to it.
    for axis in axes[:3]:
        mark_boundaries(axis)
        axis.set_xlim(1, MTDNA_LENGTH)
        axis.grid(alpha=0.15)
    axes[2].set_xlabel("NC_005089.1 position")
    fig.suptitle(f"{sample} SNV heteroplasmy QC", fontsize=15)
    fig.savefig(output_dir / f"{sample}.snv_heteroplasmy_qc.png", dpi=220)
    fig.savefig(output_dir / f"{sample}.snv_heteroplasmy_qc.pdf")
    plt.close(fig)


def plot_spectrum(pass_rows, output_dir):
    changes = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]
    samples = sorted({row["sample_id"] for row in pass_rows})
    counts = np.zeros((len(samples), len(changes)), dtype=int)
    for row in pass_rows:
        counts[samples.index(row["sample_id"]), changes.index(row["change"])] += 1
    fig, axis = plt.subplots(figsize=(12, max(3.5, 0.45 * len(samples) + 1.5)), constrained_layout=True)
    if not pass_rows:
        axis.text(0.5, 0.5, "No PASS SNVs", ha="center", va="center", transform=axis.transAxes)
        axis.set_axis_off()
        fig.savefig(output_dir / "cohort.snv_substitution_spectrum.png", dpi=220)
        fig.savefig(output_dir / "cohort.snv_substitution_spectrum.pdf")
        plt.close(fig)
        return
    left = np.zeros(len(samples))
    colors = plt.cm.tab20(np.linspace(0, 1, len(changes)))
    for index, change in enumerate(changes):
        axis.barh(samples, counts[:, index], left=left, label=change, color=colors[index])
        left += counts[:, index]
    axis.set_xlabel("PASS SNV count")
    axis.set_title("SNV substitution spectrum")
    axis.legend(frameon=False, ncol=6, fontsize=8, bbox_to_anchor=(0.5, -0.18), loc="upper center")
    fig.savefig(output_dir / "cohort.snv_substitution_spectrum.png", dpi=220)
    fig.savefig(output_dir / "cohort.snv_substitution_spectrum.pdf")
    plt.close(fig)


def plot_candidate_heatmap(pass_rows, output_dir):
    if not pass_rows:
        return
    samples = sorted({row["sample_id"] for row in pass_rows})
    variants = sorted({(int(row["standard_position"]), row["change"]) for row in pass_rows})
    matrix = np.full((len(samples), len(variants)), np.nan)
    for row in pass_rows:
        matrix[samples.index(row["sample_id"]),
               variants.index((int(row["standard_position"]), row["change"]))] = float(row["vaf_percent"])
    fig, axis = plt.subplots(figsize=(max(7, 0.55 * len(variants) + 2),
                                      max(3.5, 0.5 * len(samples) + 1.5)), constrained_layout=True)
    image = axis.imshow(matrix, aspect="auto", cmap="magma", vmin=0)
    axis.set_xticks(range(len(variants)), [f"{p}\n{c}" for p, c in variants], rotation=60, ha="right")
    axis.set_yticks(range(len(samples)), samples)
    axis.set_title("PASS SNV VAF across samples")
    fig.colorbar(image, ax=axis, label="VAF (%)")
    fig.savefig(output_dir / "cohort.pass_snv_vaf_heatmap.png", dpi=220)
    fig.savefig(output_dir / "cohort.pass_snv_vaf_heatmap.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    all_pass = []
    all_files = sorted(args.results_dir.glob("*.snv.all_sites.tsv"))
    if not all_files:
        raise FileNotFoundError(f"No *.snv.all_sites.tsv files in {args.results_dir}")
    for all_path in all_files:
        sample = all_path.name.removesuffix(".snv.all_sites.tsv")
        pass_path = args.results_dir / f"{sample}.snv.pass.tsv"
        summary_path = args.results_dir / f"{sample}.snv.summary.json"
        if not pass_path.exists() or not summary_path.exists():
            raise FileNotFoundError(f"Incomplete outputs for {sample}")
        all_rows = read_tsv(all_path)
        pass_rows = read_tsv(pass_path)
        plot_sample(sample, all_rows, pass_rows, args.output_dir)
        all_pass.extend(pass_rows)
        with summary_path.open() as handle:
            summaries.append(json.load(handle))
    plot_spectrum(all_pass, args.output_dir)
    plot_candidate_heatmap(all_pass, args.output_dir)

    summary_fields = ["sample_id", "sites_evaluated", "sites_with_any_alt",
                      "pass_snv_count", "median_depth"]
    with (args.output_dir / "cohort.snv_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)
    if all_pass:
        with (args.output_dir / "cohort.pass_snvs.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(all_pass[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(all_pass)
    print(f"Wrote cohort QC outputs to {args.output_dir}")


if __name__ == "__main__":
    main()

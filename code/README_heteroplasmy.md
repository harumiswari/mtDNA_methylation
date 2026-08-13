# Reciprocal-reference mtDNA SNV heteroplasmy pipeline

This pipeline calls SNV heteroplasmy from the same mixed-barcode reads aligned to
SphI- and NcoI-rotated mouse mtDNA references. It does not call indels.

## Analysis rule

- Use the SphI alignment as the primary evidence across most of NC_005089.1.
- Use the NcoI alignment within the configured window around the artificial
  SphI boundary (standard positions 10,758–10,759).
- Use SphI at the NcoI boundary (standard positions 9,215–9,216).
- Keep counts from the two orientations separate. Never add their depths or
  allele counts because they represent the same molecules aligned twice.
- Away from both boundaries, use the nonselected orientation as concordance QC.

## Files

1. `05_collect_dual_orientation_snv_evidence.py`
   counts A/C/G/T observations, strands, and mean base quality from both BAMs;
   performs the correct rotation separately for SphI and NcoI.
2. `06_call_snv_heteroplasmy.py`
   selects the reciprocal orientation, calculates VAF and strand bias, applies
   configurable filters, and writes all-site, candidate, PASS, and JSON summary
   outputs, plus a standards-compatible PASS VCF without a misleading diploid genotype.
3. `07_export_pass_snv_read_support.py`
   exports individual read observations for each PASS SNV from both orientations
   for IGV/read-ID auditing. Orientation rows remain separate.
4. `07_plot_snv_heteroplasmy_qc.py`
   creates per-sample depth/VAF/concordance/base-quality plots and cohort summary,
   substitution-spectrum, and PASS-VAF heatmap outputs.
5. `run_snv_heteroplasmy_pipeline.sh`
   runs the complete workflow across the configured samples.

## Requirements

```bash
python3 -m pip install pysam matplotlib numpy
samtools faidx /path/to/mt_mouse_NC_005089.1_SphI_rotated.fa
```

Both coordinate-sorted BAMs must be indexed and contain these default contigs:

```text
NC_005089.1_SphI_rotated
NC_005089.1_NcoI_rotated
```

## Run

Edit the path, filename-template, sample, and filter section at the top of
`run_snv_heteroplasmy_pipeline.sh`, then run:

```bash
bash run_snv_heteroplasmy_pipeline.sh
```

If the extracted BAM filenames differ, change only these templates:

```bash
SPHI_BAM_TEMPLATE="${INPUT_DIR}/{sample}.SphI.sorted.bam"
NCOI_BAM_TEMPLATE="${INPUT_DIR}/{sample}.NcoI.sorted.bam"
```

## Default PASS criteria

| Criterion | Default |
|---|---:|
| Selected-orientation depth | ≥100 |
| Alternate observations | ≥5 |
| VAF | 1–95% |
| Alternate support | ≥2 forward and ≥2 reverse |
| Alternate mean base quality | ≥20 |
| Fisher strand-bias p-value | ≥0.001 |
| Validation-orientation depth away from boundaries | ≥100 |
| Orientation VAF disagreement | ≤max(2 percentage points, 25% of selected VAF) |
| Reference homopolymer run | ≤6 bases |

These are starting QC thresholds, not universal biological cutoffs. Review the
VAF error floor in negative/control samples before treating 1% as a validated
limit of detection. For very deep mtDNA data, alternate-read count should also
be increased if control samples show recurrent low-frequency errors.

## Main outputs per sample

```text
SAMPLE.dual_orientation_snv_evidence.tsv
SAMPLE.snv.all_sites.tsv
SAMPLE.snv.candidates.tsv
SAMPLE.snv.pass.tsv
SAMPLE.snv.pass.vcf
SAMPLE.snv.pass_read_support.tsv
SAMPLE.snv.summary.json
```

`candidates.tsv` intentionally contains every position with at least one
alternate observation, including filtered noise. Use `pass.tsv` for the strict
candidate list and retain `all_sites.tsv` for auditability.

## Final review

Before reporting a PASS SNV, inspect both alignments in IGV and check the
read-support table for clipping, nearby indels, strand imbalance, recurrent
sequence-context errors, and read-end clustering. The pipeline annotates but
does not model all Nanopore context-specific errors, and it does not establish a
formal limit of detection without controls or orthogonal validation.

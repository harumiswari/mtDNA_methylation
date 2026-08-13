# mtDNA methylation workflow

This repository contains the ONT methylation workflow for mouse mitochondrial DNA. The pipeline starts from POD5 files, calls modified bases with Dorado, aligns the reads to dual orientated mtDNA references, extracts mtDNA-specific reads, runs modkit pileup, compares SphI and NcoI coordinate systems, and selects a final methylation signal per position.

## Overview

The workflow proceeds in four stages:

1. Modified-basecalling from POD5 data with Dorado.
2. Dual-reference alignment to SphI and NcoI composite mtDNA references.
3. mtDNA extraction and modkit pileup in each orientation.
4. Coordinate conversion, reciprocal comparison, and final orientation selection.

## Why two orientations?

The same mtDNA molecule is represented in two rotated coordinate systems:

- SphI-rotated mtDNA reference
- NcoI-rotated mtDNA reference

These references are intentionally separated because the same molecule can be represented differently depending on the linearization point. Keeping the alignments separate allows the workflow to avoid coordinate ambiguity and to evaluate orientation-specific methylation behavior near restriction boundaries.

## Repository structure

```text
code/
  01_dna_basemod_call_batch.sh
  02_modkit_pileup_batch.sh
  03_compare_modkit_orientations.py
  04_select_final_modkit_orientation(2).py
Ref/
  GRCm39_nuclear_SphI_mtDNA_DCS.fa
  GRCm39_nuclear_NcoI_mtDNA_DCS.fa
  mt_mouse_NC_005089.1_SphI_rotated.fa
  mt_mouse_NC_005089.1_NcoI_rotated.fa
```

## Requirements

Activate the analysis environment:

```bash
conda activate longread
```

Required software:

- Dorado 2.0.1
- samtools
- modkit
- Python 3
- pysam (for downstream processing)

Expected Dorado paths:

```text
/home/harumiswari/dorado-2.0.1-linux-x64/bin/dorado
/home/harumiswari/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```

## Step 1: modified basecalling and dual-reference alignment

Script: `code/01_dna_basemod_call_batch.sh`

This script does the following for each sample (`XX_01` through `XX_08`):

1. Uses Dorado basecaller with modified-base support for `5mCG_5hmCG`.
2. Aligns the basecalled reads to the SphI composite reference.
3. Aligns the same reads to the NcoI composite reference.
4. Extracts primary reads assigned to the mtDNA contig in each alignment.
5. Filters reads by:
   - MAPQ >= 20
   - excluding unmapped, secondary, supplementary, duplicate, QCFail reads
6. Sorts and indexes the retained mtDNA BAMs.

Output directory:

```text
XX_methylation_basemod_results/
```

Representative outputs:

```text
XX_01.basecalled_mods.bam
XX_01.SphI.composite.sorted.bam
XX_01.NcoI.composite.sorted.bam
XX_01.SphI.sorted.bam
XX_01.NcoI.sorted.bam
```

## Step 2: modkit pileup in each orientation

Script: `code/02_modkit_pileup_batch.sh`

This script runs modkit pileup separately on the SphI and NcoI mtDNA BAMs against their matching mtDNA-only references:

```text
SphI BAM -> mt_mouse_NC_005089.1_SphI_rotated.fa
NcoI BAM -> mt_mouse_NC_005089.1_NcoI_rotated.fa
```

It requests:

- 5mC
- 5hmC
- CpG context
- combined strands
- BGZF-compressed output

Output directory:

```text
XX_methylation_modkit_results/
```

Example outputs:

```text
XX_01.SphI.5mC_5hmC.CpG.bed.gz
XX_01.NcoI.5mC_5hmC.CpG.bed.gz
```

## Step 3: convert coordinates and compare SphI/NcoI methylation

Script: `code/03_compare_modkit_orientations.py`

This script:

1. Reads the SphI and NcoI modkit BED.gz outputs.
2. Converts rotated positions back to standard NC_005089.1 coordinates.
3. Preserves 5mC and 5hmC counts and fractions.
4. Calculates total modified fraction from counts and coverage.
5. Flags positions near the restriction boundaries.
6. Writes a per-sample comparison table.
7. Writes a multi-sample summary table.

Key logic:

- SphI reference uses `SPHI_FIRST_STANDARD_POSITION = 10_759`
- NcoI conversion is handled through the reciprocal coordinate mapping
- positions near the artificial boundaries are flagged as `boundary_review`

Output directory:

```text
XX_methylation_comparison_results/
```

Representative output:

```text
XX_01.SphI_NcoI.comparison.tsv
XX_replicate_summary.tsv
```

## Step 4: choose the final methylation orientation

Script: `code/04_select_final_modkit_orientation(2).py`

This script selects the final per-position methylation call using the comparison table. The decision logic is:

- if a site is near the SphI boundary, prefer NcoI when coverage is available
- if a site is near the NcoI boundary, prefer SphI when coverage is available
- if one orientation is missing, rescue with the other orientation
- otherwise prefer the orientation with higher coverage

It also applies a QC gate based on:

- minimum coverage threshold
- minimum relative coverage
- maximum orientation disagreement threshold

The script writes:

- a final per-position methylation table per sample
- a summary table across samples

Output directory:

```text
XX_methylation_final_results/
```

Representative outputs:

```text
XX_01.final_methylation.tsv
XX_final_methylation_summary.tsv
```

## Boundary interpretation

The restriction-site neighborhoods are critical for interpretation:

- SphI boundary: standard positions 10,758-10,759
- NcoI boundary: standard position 9,216
- default boundary window: +/-100 bp

Positions near these boundaries are flagged for special review because the read population can be biased by the artificial linearization point and by an end-related coverage effect.

## QC recommendations

Before interpreting methylation calls, review:

- mtDNA coverage along the full 16,299 bp molecule
- read assignment quality to mtDNA contigs
- SphI/NcoI agreement away from the artificial boundaries
- boundary-localized disagreement
- read-end and clipping patterns
- shifts in percentage estimates when switching orientation

## Notes

This workflow is designed to be upstream of the reciprocal-orientation heteroplasmy analysis. It generates the cleaned mtDNA BAMs and orientation-specific methylation tables required for downstream heteroplasmy comparisons and validation.

The methylation analysis is intentionally separated from the heteroplasmy analysis: the first workflow handles modified-basecalling and methylation assessment, while the heteroplasmy repository handles reciprocal BAM evidence and SNV-level variant calling.

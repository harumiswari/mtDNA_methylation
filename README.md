# Mouse mtDNA ONT Methylation Workflow

This workflow analyzes Oxford Nanopore POD5 data for mouse mitochondrial DNA (mtDNA) 5mC and 5hmC modification patterns across eight future biological replicates.

## Workflow Overview

```text
POD5 files
    |
    v
Dorado modified-basecalling
    |
    v
Modification-tagged BAM
    |
    +--> GRCm39 nuclear + SphI-mtDNA + DCS composite alignment
    |        |
    |        +--> primary, MAPQ-filtered SphI mtDNA BAM
    |
    +--> GRCm39 nuclear + NcoI-mtDNA + DCS composite alignment
             |
             +--> primary, MAPQ-filtered NcoI mtDNA BAM
                          |
                          v
                    modkit pileup
                          |
                          v
             SphI/NcoI coordinate comparison
```

The two composite alignments are intentional. The same mtDNA sequence is represented in SphI and NcoI orientations, so both orientations must not be placed in one reference. Separate alignments allow NUMT sequences to compete against mtDNA while keeping the coordinate systems unambiguous.

## Requirements

Activate the environment containing the required tools:

```bash
conda activate longread
```

Required software:

- Dorado 2.0.1
- samtools
- modkit
- Python 3

The workflow expects the following Dorado paths:

```text
Dorado: /home/harumiswari/dorado-2.0.1-linux-x64/bin/dorado
Model:  /home/harumiswari/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```

## References

References are stored in `Ref/`:

```text
Ref/GRCm39_nuclear_SphI_mtDNA_DCS.fa
Ref/GRCm39_nuclear_NcoI_mtDNA_DCS.fa
Ref/mt_mouse_NC_005089.1_SphI_rotated.fa
Ref/mt_mouse_NC_005089.1_NcoI_rotated.fa
```

The composite references contain:

- GRCm39 nuclear contigs without chrM
- One orientation-specific mtDNA contig
- The DCS lambda control contig

The SphI and NcoI composite references must remain separate. The mtDNA-only references are used by `modkit` after the composite alignment has assigned reads to mtDNA.

## Sample Layout

The batch scripts use placeholder sample IDs `XX_01` through `XX_08`.

Create one POD5 directory per sample:

```text
/home/harumiswari/Documents/mtdna/pod5_XX_01
/home/harumiswari/Documents/mtdna/pod5_XX_02
/home/harumiswari/Documents/mtdna/pod5_XX_03
/home/harumiswari/Documents/mtdna/pod5_XX_04
/home/harumiswari/Documents/mtdna/pod5_XX_05
/home/harumiswari/Documents/mtdna/pod5_XX_06
/home/harumiswari/Documents/mtdna/pod5_XX_07
/home/harumiswari/Documents/mtdna/pod5_XX_08
```

Each directory should contain one or more `.pod5` files.

## Step 1: Modified-Basecalling and Alignment

Run:

```bash
bash code/01_dna_basemod_call_batch.sh
```

For each sample, the script:

1. Basecalls POD5 once using Dorado SUP.
2. Detects `5mCG` and `5hmCG`.
3. Aligns the tagged reads to the SphI composite reference.
4. Aligns the same tagged reads to the NcoI composite reference.
5. Extracts the matching mtDNA contig from each composite BAM.
6. Excludes unmapped, secondary, supplementary, QC-failed, and duplicate reads.
7. Requires MAPQ >= 20 for extracted mtDNA reads.
8. Sorts, checks, and indexes the BAM files.

Main output directory:

```text
XX_methylation_basemod_results/
```

Important outputs per sample:

```text
XX_01.SphI.composite.sorted.bam
XX_01.NcoI.composite.sorted.bam
XX_01.SphI.sorted.bam
XX_01.NcoI.sorted.bam
XX_01.basecalled_mods.bam
```

The composite BAMs are retained for NUMT and alignment-quality review. The mtDNA-only BAMs are used for `modkit`.

## Step 2: modkit Pileup

Run:

```bash
bash code/02_modkit_pileup_batch.sh
```

The script runs `modkit pileup` separately on the SphI and NcoI mtDNA BAMs using the matching mtDNA-only reference:

```text
SphI BAM -> mt_mouse_NC_005089.1_SphI_rotated.fa
NcoI BAM -> mt_mouse_NC_005089.1_NcoI_rotated.fa
```

It requests:

- 5mC
- 5hmC
- CpG sites
- combined strands
- BGZF-compressed output

Output directory:

```text
XX_methylation_modkit_results/
```

Expected files per sample:

```text
XX_01.SphI.5mC_5hmC.CpG.bed.gz
XX_01.NcoI.5mC_5hmC.CpG.bed.gz
```

## Step 3: Coordinate Conversion and Comparison

Run:

```bash
python code/03_compare_modkit_orientations.py
```

This script:

1. Reads the SphI and NcoI modkit outputs.
2. Converts both rotated coordinate systems to standard NC_005089.1 coordinates.
3. Preserves 5mC and 5hmC counts and fractions.
4. Calculates total modified fraction from counts and coverage.
5. Flags positions near the SphI and NcoI linearization boundaries.
6. Writes one comparison table per sample.
7. Writes an eight-replicate summary table.

Output directory:

```text
XX_methylation_comparison_results/
```

The combined methylation fraction is count-weighted within each orientation. Percentages from SphI and NcoI are not averaged directly.

## Boundary Interpretation

The restriction-site neighborhoods require special review:

- SphI boundary: approximately standard positions 10,758-10,759.
- NcoI boundary: approximately standard position 9,216.
- The comparison script flags a default window of +/-100 bp.

Around a restriction site, use the opposite enzyme orientation as the preferred evidence because molecules cut by that enzyme may show end-related coverage or modification-call bias near the cut.

## Quality Control

Before interpreting methylation results, check:

```bash
samtools quickcheck XX_methylation_basemod_results/XX_01.SphI.sorted.bam
samtools quickcheck XX_methylation_basemod_results/XX_01.NcoI.sorted.bam
samtools idxstats XX_methylation_basemod_results/XX_01.SphI.sorted.bam
samtools idxstats XX_methylation_basemod_results/XX_01.NcoI.sorted.bam
```

Review:

- mtDNA coverage across the full 16,299 bp molecule
- MAPQ distribution
- aligned read length and read-end concentration
- supplementary alignment and SA-tag patterns in the composite BAMs
- agreement between SphI and NcoI orientations
- replicate-level consistency

Do not merge the SphI and NcoI BAMs. They use different coordinate systems and represent the same molecules in alternative reference orientations.

## Current Status

The scripts are prepared for future samples named `XX_01` through `XX_08`. No sample analysis will run until the corresponding POD5 directories and files exist.

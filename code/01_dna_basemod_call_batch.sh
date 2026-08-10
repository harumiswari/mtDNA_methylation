#!/usr/bin/env bash
set -euo pipefail

SAMPLES=(XX_{01..08})

DORADO="/home/harumiswari/dorado-2.0.1-linux-x64/bin/dorado"
MODEL="/home/harumiswari/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

SPHI_REF="/home/harumiswari/Documents/mtdna/Ref/GRCm39_nuclear_SphI_mtDNA_DCS.fa"
NCOI_REF="/home/harumiswari/Documents/mtdna/Ref/GRCm39_nuclear_NcoI_mtDNA_DCS.fa"

SPHI_CONTIG="NC_005089.1_SphI_rotated"
NCOI_CONTIG="NC_005089.1_NcoI_rotated"

POD5_ROOT="/home/harumiswari/Documents/mtdna"
OUTPUT_ROOT="/home/harumiswari/Documents/mtdna/XX_methylation_basemod_results"

THREADS=8
MIN_MAPQ=20
EXCLUDE_FLAGS=2308

mkdir -p "$OUTPUT_ROOT"

for SAMPLE in "${SAMPLES[@]}"; do

    POD5_DIR="${POD5_ROOT}/pod5_${SAMPLE}"
    BASECALL_BAM="${OUTPUT_ROOT}/${SAMPLE}.basecalled_mods.bam"

    SPHI_COMPOSITE_BAM="${OUTPUT_ROOT}/${SAMPLE}.SphI.composite.sorted.bam"
    NCOI_COMPOSITE_BAM="${OUTPUT_ROOT}/${SAMPLE}.NcoI.composite.sorted.bam"

    SPHI_BAM="${OUTPUT_ROOT}/${SAMPLE}.SphI.sorted.bam"
    NCOI_BAM="${OUTPUT_ROOT}/${SAMPLE}.NcoI.sorted.bam"

    echo "Processing ${SAMPLE}"

    # Basecall once with 5mC and 5hmC detection
    "$DORADO" basecaller \
        "$MODEL" \
        "$POD5_DIR" \
        --modified-bases 5mCG_5hmCG \
        > "$BASECALL_BAM" \
        2> "${OUTPUT_ROOT}/${SAMPLE}.basecaller.log"

    # Align to the SphI composite reference.
    "$DORADO" aligner \
        "$SPHI_REF" \
        "$BASECALL_BAM" \
        --threads "$THREADS" \
        2> "${OUTPUT_ROOT}/${SAMPLE}.SphI.aligner.log" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$SPHI_COMPOSITE_BAM"

    samtools index -@ "$THREADS" "$SPHI_COMPOSITE_BAM"

    # Extract only reads assigned to the SphI mtDNA contig.
    samtools view \
        -@ "$THREADS" \
        -b \
        -F "$EXCLUDE_FLAGS" \
        -q "$MIN_MAPQ" \
        -o "$SPHI_BAM" \
        "$SPHI_COMPOSITE_BAM" \
        "$SPHI_CONTIG"

    samtools quickcheck -v "$SPHI_BAM"
    samtools index -@ "$THREADS" "$SPHI_BAM"

    # Align to the NcoI composite reference.
    "$DORADO" aligner \
        "$NCOI_REF" \
        "$BASECALL_BAM" \
        --threads "$THREADS" \
        2> "${OUTPUT_ROOT}/${SAMPLE}.NcoI.aligner.log" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$NCOI_COMPOSITE_BAM"

    samtools index -@ "$THREADS" "$NCOI_COMPOSITE_BAM"

    # Extract primary, confidently assigned NcoI mtDNA reads.
    samtools view \
        -@ "$THREADS" \
        -b \
        -F "$EXCLUDE_FLAGS" \
        -q "$MIN_MAPQ" \
        -o "$NCOI_BAM" \
        "$NCOI_COMPOSITE_BAM" \
        "$NCOI_CONTIG"

    samtools quickcheck -v "$NCOI_BAM"
    samtools index -@ "$THREADS" "$NCOI_BAM"

    echo "Completed ${SAMPLE}"

done

echo "All samples completed."

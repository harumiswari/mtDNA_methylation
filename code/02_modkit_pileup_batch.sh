#!/usr/bin/env bash
set -euo pipefail

SAMPLES=(XX_{01..08})

REFERENCE_DIR="/home/harumiswari/Documents/mtdna/Ref"
SPHI_REF="${REFERENCE_DIR}/mt_mouse_NC_005089.1_SphI_rotated.fa"
NCOI_REF="${REFERENCE_DIR}/mt_mouse_NC_005089.1_NcoI_rotated.fa"

INPUT_DIR="/home/harumiswari/Documents/mtdna/XX_methylation_basemod_results"
OUTPUT_DIR="/home/harumiswari/Documents/mtdna/XX_methylation_modkit_results"

mkdir -p "$OUTPUT_DIR"

for SAMPLE in "${SAMPLES[@]}"; do

    SPHI_BAM="${INPUT_DIR}/${SAMPLE}.SphI.sorted.bam"
    NCOI_BAM="${INPUT_DIR}/${SAMPLE}.NcoI.sorted.bam"

    SPHI_OUTPUT="${OUTPUT_DIR}/${SAMPLE}.SphI.5mC_5hmC.CpG.bed.gz"
    NCOI_OUTPUT="${OUTPUT_DIR}/${SAMPLE}.NcoI.5mC_5hmC.CpG.bed.gz"

    echo "Processing ${SAMPLE}"

    # Call 5mC and 5hmC on the SphI-oriented alignment.
    modkit pileup \
        "$SPHI_BAM" \
        "$SPHI_OUTPUT" \
        --modified-bases 5mC 5hmC \
        --cpg \
        --combine-strands \
        --ref "$SPHI_REF" \
        --bgzf \
        --log "${OUTPUT_DIR}/${SAMPLE}.SphI.modkit.log"

    # Call 5mC and 5hmC on the NcoI-oriented alignment.
    modkit pileup \
        "$NCOI_BAM" \
        "$NCOI_OUTPUT" \
        --modified-bases 5mC 5hmC \
        --cpg \
        --combine-strands \
        --ref "$NCOI_REF" \
        --bgzf \
        --log "${OUTPUT_DIR}/${SAMPLE}.NcoI.modkit.log"

    echo "Completed ${SAMPLE}"
done

echo "All samples completed."
echo "Results: $OUTPUT_DIR"

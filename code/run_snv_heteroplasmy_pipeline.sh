#!/usr/bin/env bash
set -euo pipefail

# Edit these paths once. BAM templates must contain the literal token {sample}.
SAMPLES=(XX_{01..08})
INPUT_DIR="/home/harumiswari/Documents/mtdna/XX_methylation_basemod_results"
RESULTS_DIR="/home/harumiswari/Documents/mtdna/XX_snv_heteroplasmy_results"
SPHI_REFERENCE="/home/harumiswari/Documents/mtdna/Ref/mt_mouse_NC_005089.1_SphI_rotated.fa"
SPHI_BAM_TEMPLATE="${INPUT_DIR}/{sample}.SphI.sorted.bam"
NCOI_BAM_TEMPLATE="${INPUT_DIR}/{sample}.NcoI.sorted.bam"

# Read/base filters used while collecting evidence.
MIN_MAPQ=20
MIN_BASEQ=20
BOUNDARY_WINDOW=100

# Candidate filters. VAF values are fractions: 0.01 = 1%.
MIN_DEPTH=100
MIN_ALT_COUNT=5
MIN_VAF=0.01
MAX_VAF=0.95
MIN_ALT_PER_STRAND=2
MIN_ALT_MEAN_BASEQ=20
STRAND_BIAS_P=0.001
MIN_VALIDATION_DEPTH=100
MAX_ORIENTATION_VAF_DIFF=0.02
MAX_ORIENTATION_RELATIVE_DIFF=0.25
MAX_HOMOPOLYMER_RUN=6

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COLLECTOR="${SCRIPT_DIR}/05_collect_dual_orientation_snv_evidence.py"
CALLER="${SCRIPT_DIR}/06_call_snv_heteroplasmy.py"
PLOTTER="${SCRIPT_DIR}/07_plot_snv_heteroplasmy_qc.py"
SUPPORT_EXPORTER="${SCRIPT_DIR}/07_export_pass_snv_read_support.py"

mkdir -p "$RESULTS_DIR"
for SAMPLE in "${SAMPLES[@]}"; do
    SPHI_BAM="${SPHI_BAM_TEMPLATE//\{sample\}/$SAMPLE}"
    NCOI_BAM="${NCOI_BAM_TEMPLATE//\{sample\}/$SAMPLE}"
    EVIDENCE="${RESULTS_DIR}/${SAMPLE}.dual_orientation_snv_evidence.tsv"
    PREFIX="${RESULTS_DIR}/${SAMPLE}.snv"

    echo "Collecting SNV evidence for ${SAMPLE}"
    python3 "$COLLECTOR" \
        --sample "$SAMPLE" \
        --sphi-bam "$SPHI_BAM" \
        --ncoi-bam "$NCOI_BAM" \
        --sphi-reference "$SPHI_REFERENCE" \
        --output "$EVIDENCE" \
        --min-mapq "$MIN_MAPQ" \
        --min-baseq "$MIN_BASEQ" \
        --boundary-window "$BOUNDARY_WINDOW"

    echo "Calling SNV heteroplasmy for ${SAMPLE}"
    python3 "$CALLER" \
        --evidence "$EVIDENCE" \
        --output-prefix "$PREFIX" \
        --min-depth "$MIN_DEPTH" \
        --min-alt-count "$MIN_ALT_COUNT" \
        --min-vaf "$MIN_VAF" \
        --max-vaf "$MAX_VAF" \
        --min-alt-per-strand "$MIN_ALT_PER_STRAND" \
        --min-alt-mean-baseq "$MIN_ALT_MEAN_BASEQ" \
        --strand-bias-p "$STRAND_BIAS_P" \
        --min-validation-depth "$MIN_VALIDATION_DEPTH" \
        --max-orientation-vaf-diff "$MAX_ORIENTATION_VAF_DIFF" \
        --max-orientation-relative-diff "$MAX_ORIENTATION_RELATIVE_DIFF" \
        --max-homopolymer-run "$MAX_HOMOPOLYMER_RUN"

    python3 "$SUPPORT_EXPORTER" \
        --pass-calls "${PREFIX}.pass.tsv" \
        --sphi-bam "$SPHI_BAM" \
        --ncoi-bam "$NCOI_BAM" \
        --output "${PREFIX}.pass_read_support.tsv" \
        --min-mapq "$MIN_MAPQ" \
        --min-baseq "$MIN_BASEQ"
done

python3 "$PLOTTER" \
    --results-dir "$RESULTS_DIR" \
    --output-dir "${RESULTS_DIR}/plots"

echo "SNV heteroplasmy pipeline complete: ${RESULTS_DIR}"

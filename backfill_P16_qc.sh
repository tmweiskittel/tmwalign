#!/usr/bin/env bash

set -euo pipefail


# ============================================================================
# Configuration
# ============================================================================

SAMPLE="P16"

REPO_PATH="/home/jupyter/repos/tmw_align_emseq"

WORK_DIR="/home/jupyter/data/p16_rebuild"

FASTP_JSON_DIR="/home/jupyter/data/qc/fastp"

COVERAGE_SCRIPT="${REPO_PATH}/workflow/scripts/coverage_qc.py"
SUMMARY_SCRIPT="${REPO_PATH}/workflow/scripts/sample_qc_summary_once.py"

RAW_R1_GS="gs://weiskittel-projects1/radnecrosis/raw_data/01.RawData/P16/P16_WKDL250008855-1A_2353KGLT4_L5_1.fq.gz"
RAW_R2_GS="gs://weiskittel-projects1/radnecrosis/raw_data/01.RawData/P16/P16_WKDL250008855-1A_2353KGLT4_L5_2.fq.gz"

PROCESSED_ROOT="gs://weiskittel-projects1/radnecrosis/processed/P16"

GENOME_SIZE=3100000000
MIN_MAPQ=30
MIN_BASEQ=0
EXCLUDED_CONTIGS="lambda,pUC19"


# ============================================================================
# Paths
# ============================================================================

mkdir -p "${WORK_DIR}"
mkdir -p "${FASTP_JSON_DIR}"

RAW_R1="${WORK_DIR}/P16.R1.fastq.gz"
RAW_R2="${WORK_DIR}/P16.R2.fastq.gz"

TRIMMED_R1="${WORK_DIR}/P16.trimmed.R1.fastq.gz"
TRIMMED_R2="${WORK_DIR}/P16.trimmed.R2.fastq.gz"

FASTP_JSON="${FASTP_JSON_DIR}/P16.fastp.json"
FASTP_HTML="${WORK_DIR}/P16.fastp.html"

RAW_BAM="${WORK_DIR}/P16.aligned.sorted.bam"

FINAL_BAM="${WORK_DIR}/P16.aligned.sorted.filt.bl.bam"
FINAL_BAI="${WORK_DIR}/P16.aligned.sorted.filt.bl.bam.bai"

CPG_FILE="${WORK_DIR}/P16.CpG.methylKit.gz"
LAMBDA_QC="${WORK_DIR}/P16.lambda_qc.tsv"

COVERAGE_QC="${WORK_DIR}/P16.coverage_qc.tsv"
SUMMARY_QC="${WORK_DIR}/P16.qc_summary.tsv"

FASTP_LOG="${WORK_DIR}/P16.fastp.log"
COVERAGE_LOG="${WORK_DIR}/P16.coverage_qc.log"
SUMMARY_LOG="${WORK_DIR}/P16.qc_summary.log"


# ============================================================================
# Validate required programs and scripts
# ============================================================================

for program in gcloud fastp python3 samtools; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: ${program}" >&2
        exit 1
    fi
done

if [[ ! -f "${COVERAGE_SCRIPT}" ]]; then
    echo "ERROR: Missing coverage script:"
    echo "  ${COVERAGE_SCRIPT}"
    exit 1
fi

if [[ ! -f "${SUMMARY_SCRIPT}" ]]; then
    echo "ERROR: Missing summary script:"
    echo "  ${SUMMARY_SCRIPT}"
    exit 1
fi


# ============================================================================
# Step 1: Download raw FASTQs
# ============================================================================

echo "Downloading raw FASTQs..."

gcloud storage cp \
    "${RAW_R1_GS}" \
    "${RAW_R1}"

gcloud storage cp \
    "${RAW_R2_GS}" \
    "${RAW_R2}"


# ============================================================================
# Step 2: Recreate fastp metrics exactly as the original pipeline did
# ============================================================================

echo "Running fastp..."

fastp \
    --in1 "${RAW_R1}" \
    --in2 "${RAW_R2}" \
    --out1 "${TRIMMED_R1}" \
    --out2 "${TRIMMED_R2}" \
    --html "${FASTP_HTML}" \
    --json "${FASTP_JSON}" \
    --thread 4 \
    > "${FASTP_LOG}" 2>&1

if [[ ! -s "${FASTP_JSON}" ]]; then
    echo "ERROR: fastp JSON was not created." >&2
    exit 1
fi

echo "fastp JSON created:"
echo "  ${FASTP_JSON}"


# Raw and trimmed FASTQs are no longer needed.
rm -f \
    "${RAW_R1}" \
    "${RAW_R2}" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}"


# ============================================================================
# Step 3: Download preserved analysis artifacts
# ============================================================================

echo "Downloading preserved analysis artifacts..."

gcloud storage cp \
    "${PROCESSED_ROOT}/P16.aligned.sorted.bam" \
    "${RAW_BAM}"

gcloud storage cp \
    "${PROCESSED_ROOT}/P16.aligned.sorted.filt.bl.bam" \
    "${FINAL_BAM}"

gcloud storage cp \
    "${PROCESSED_ROOT}/P16.aligned.sorted.filt.bl.bam.bai" \
    "${FINAL_BAI}"

gcloud storage cp \
    "${PROCESSED_ROOT}/P16.CpG.methylKit.gz" \
    "${CPG_FILE}"

gcloud storage cp \
    "${PROCESSED_ROOT}/P16.lambda_qc.tsv" \
    "${LAMBDA_QC}"


# ============================================================================
# Step 4: Validate BAMs
# ============================================================================

echo "Validating BAMs..."

samtools quickcheck \
    "${RAW_BAM}" \
    "${FINAL_BAM}"

samtools idxstats \
    "${FINAL_BAM}" \
    >/dev/null


# ============================================================================
# Step 5: Generate new coverage QC
# ============================================================================

echo "Generating corrected coverage QC..."

python3 "${COVERAGE_SCRIPT}" \
    --sample "${SAMPLE}" \
    --cpg "${CPG_FILE}" \
    --bam "${FINAL_BAM}" \
    --fastp-json "${FASTP_JSON}" \
    --output "${COVERAGE_QC}" \
    --genome-size "${GENOME_SIZE}" \
    --min-mapq "${MIN_MAPQ}" \
    --min-baseq "${MIN_BASEQ}" \
    --excluded-contigs "${EXCLUDED_CONTIGS}" \
    > "${COVERAGE_LOG}" 2>&1

if [[ ! -s "${COVERAGE_QC}" ]]; then
    echo "ERROR: Coverage QC was not generated." >&2
    exit 1
fi


# ============================================================================
# Step 6: Rebuild QC summary
# ============================================================================

echo "Rebuilding QC summary..."

python3 "${SUMMARY_SCRIPT}" \
    --sample "${SAMPLE}" \
    --fastp-json "${FASTP_JSON}" \
    --raw-bam "${RAW_BAM}" \
    --final-bam "${FINAL_BAM}" \
    --lambda-qc "${LAMBDA_QC}" \
    --coverage-qc "${COVERAGE_QC}" \
    --output "${SUMMARY_QC}" \
    > "${SUMMARY_LOG}" 2>&1

if [[ ! -s "${SUMMARY_QC}" ]]; then
    echo "ERROR: QC summary was not generated." >&2
    exit 1
fi


# ============================================================================
# Step 7: Check for obsolete schema fields
# ============================================================================

echo "Validating new QC summary schema..."

if head -n 1 "${SUMMARY_QC}" \
    | grep -qE \
        'estimated_genome_depth|cpg_sites_ge_10x|mean_coverage($|\t)'; then

    echo "ERROR: New QC summary still contains obsolete fields." >&2
    exit 1
fi

for required_field in \
    raw_fastq_coverage \
    mean_aligned_base_coverage \
    mean_coverage_called_cpgs \
    median_coverage_called_cpgs \
    coverage_weighted_methylation_fraction
do
    if ! head -n 1 "${SUMMARY_QC}" \
        | tr '\t' '\n' \
        | grep -qx "${required_field}"; then

        echo "ERROR: Missing required field: ${required_field}" >&2
        exit 1
    fi
done


# ============================================================================
# Step 8: Upload corrected QC into native processed/P16 location
# ============================================================================

echo "Uploading corrected native QC outputs..."

gcloud storage cp \
    "${COVERAGE_QC}" \
    "${PROCESSED_ROOT}/P16.coverage_qc.tsv"

gcloud storage cp \
    "${SUMMARY_QC}" \
    "${PROCESSED_ROOT}/P16.qc_summary.tsv"


# ============================================================================
# Step 9: Verify uploaded native summary
# ============================================================================

echo "Verifying uploaded QC summary..."

uploaded_header=$(
    gcloud storage cat \
        "${PROCESSED_ROOT}/P16.qc_summary.tsv" \
        | head -n 1
)

if printf '%s\n' "${uploaded_header}" \
    | grep -qE \
        'estimated_genome_depth|cpg_sites_ge_10x|mean_coverage($|\t)'; then

    echo "ERROR: Uploaded P16 summary appears to have the old schema." >&2
    exit 1
fi


# ============================================================================
# Step 10: Remove stale analysis-side P16 cache
# ============================================================================

ANALYSIS_CACHE="/home/jupyter/repos/tmw_analysis_emseq/qc/P16.qc_summary.tsv"

if [[ -f "${ANALYSIS_CACHE}" ]]; then
    echo "Removing stale analysis-side P16 QC cache..."
    rm -f "${ANALYSIS_CACHE}"
fi


# ============================================================================
# Step 11: Remove large temporary BAM files
# ============================================================================

echo "Removing temporary BAM and methylation inputs..."

rm -f \
    "${RAW_BAM}" \
    "${FINAL_BAM}" \
    "${FINAL_BAI}" \
    "${CPG_FILE}" \
    "${LAMBDA_QC}"


# ============================================================================
# Complete
# ============================================================================

echo
echo "P16 QC backfill complete."
echo
echo "Retained local artifacts:"
echo "  ${FASTP_JSON}"
echo "  ${FASTP_HTML}"
echo "  ${COVERAGE_QC}"
echo "  ${SUMMARY_QC}"
echo "  ${FASTP_LOG}"
echo "  ${COVERAGE_LOG}"
echo "  ${SUMMARY_LOG}"
echo
echo "Native cloud outputs:"
echo "  ${PROCESSED_ROOT}/P16.coverage_qc.tsv"
echo "  ${PROCESSED_ROOT}/P16.qc_summary.tsv"

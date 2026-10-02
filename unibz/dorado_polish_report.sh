#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Dorado polishing report
#
# Usage:
#
#   ./dorado_polish_report.sh \
#       --draft consensus_assembly.fasta \
#       --polished consensus.fastq \
#       --vcf variants.vcf \
#       --output polishing_report.txt
#
# Optional:
#
#   --quast-dir quast_polishing_comparison
#
# Requirements:
#   Required:
#       seqkit OR seqtk   (only needed if polished input is FASTQ)
#       bcftools           (for VCF statistics)
#
#   Optional:
#       dnadiff            (from MUMmer)
#       quast.py           (for assembly statistics)
#
###############################################################################

usage() {
    cat <<EOF

Usage:
  $0 --draft DRAFT.fa --polished POLISHED.fa/fq --vcf VARIANTS.vcf --output REPORT.txt

Options:
  --draft       Original draft assembly FASTA
  --polished    Dorado polished assembly FASTA or FASTQ
  --vcf         Dorado variants.vcf
  --output      Output report filename

  --quast-dir   Optional directory for QUAST results
  --help        Show this help

Example:
  $0 \\
      --draft consensus_assembly.fasta \\
      --polished consensus.fastq \\
      --vcf variants.vcf \\
      --output dorado_polishing_report.txt

EOF
}

DRAFT=""
POLISHED=""
VCF=""
REPORT=""
QUAST_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --draft)
            DRAFT="$2"
            shift 2
            ;;
        --polished)
            POLISHED="$2"
            shift 2
            ;;
        --vcf)
            VCF="$2"
            shift 2
            ;;
        --output)
            REPORT="$2"
            shift 2
            ;;
        --quast-dir)
            QUAST_DIR="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

###############################################################################
# Check arguments
###############################################################################

if [[ -z "$DRAFT" || -z "$POLISHED" || -z "$VCF" || -z "$REPORT" ]]; then
    echo "ERROR: --draft, --polished, --vcf and --output are required."
    usage
    exit 1
fi

for file in "$DRAFT" "$POLISHED" "$VCF"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done

###############################################################################
# Temporary directory
###############################################################################

TMPDIR_REPORT=$(mktemp -d)

cleanup() {
    rm -rf "$TMPDIR_REPORT"
}

trap cleanup EXIT

###############################################################################
# Convert polished FASTQ -> FASTA if necessary
###############################################################################

POLISHED_FASTA="$POLISHED"

case "$POLISHED" in
    *.fastq|*.fq|*.fastq.gz|*.fq.gz)

        POLISHED_FASTA="$TMPDIR_REPORT/polished.fasta"

        if command -v seqkit >/dev/null 2>&1; then

            echo "Converting polished FASTQ to FASTA using seqkit..."

            seqkit fq2fa "$POLISHED" \
                -o "$POLISHED_FASTA"

        elif command -v seqtk >/dev/null 2>&1; then

            echo "Converting polished FASTQ to FASTA using seqtk..."

            if [[ "$POLISHED" == *.gz ]]; then
                zcat "$POLISHED" | seqtk seq -a - > "$POLISHED_FASTA"
            else
                seqtk seq -a "$POLISHED" > "$POLISHED_FASTA"
            fi

        else
            echo "ERROR: polished assembly is FASTQ but neither seqkit nor seqtk is available."
            exit 1
        fi
        ;;
esac

###############################################################################
# Helper functions
###############################################################################

fasta_length() {
    awk '
        BEGIN { n=0 }
        /^>/ { next }
        { gsub(/[ \t\r\n]/,""); n += length($0) }
        END { print n }
    ' "$1"
}

fasta_contigs() {
    grep -c '^>' "$1"
}

fasta_gc() {
    awk '
        BEGIN { gc=0; total=0 }
        /^>/ { next }
        {
            gsub(/[ \t\r\n]/,"")
            seq=toupper($0)

            gc += gsub(/[GC]/,"",seq)
            total += length($0)
        }
        END {
            if (total > 0)
                printf "%.3f", (gc / total) * 100
            else
                print "NA"
        }
    ' "$1"
}

###############################################################################
# Assembly statistics
###############################################################################

DRAFT_LENGTH=$(fasta_length "$DRAFT")
POLISHED_LENGTH=$(fasta_length "$POLISHED_FASTA")

DRAFT_CONTIGS=$(fasta_contigs "$DRAFT")
POLISHED_CONTIGS=$(fasta_contigs "$POLISHED_FASTA")

DRAFT_GC=$(fasta_gc "$DRAFT")
POLISHED_GC=$(fasta_gc "$POLISHED_FASTA")

LENGTH_DIFF=$((POLISHED_LENGTH - DRAFT_LENGTH))

###############################################################################
# VCF statistics
###############################################################################

VCF_TOTAL="NA"
VCF_SNP="NA"
VCF_INS="NA"
VCF_DEL="NA"
VCF_OTHER="NA"

if command -v bcftools >/dev/null 2>&1; then

    VCF_TOTAL=$(bcftools view -H "$VCF" | wc -l)

    VCF_SNP=$(bcftools view -v snps -H "$VCF" | wc -l)

    VCF_INS=$(bcftools view -v indels -H "$VCF" \
        | awk '
            {
                ref=$4
                alt=$5

                if (length(alt) > length(ref))
                    n++
            }
            END { print n+0 }
        ')

    VCF_DEL=$(bcftools view -v indels -H "$VCF" \
        | awk '
            {
                ref=$4
                alt=$5

                if (length(alt) < length(ref))
                    n++
            }
            END { print n+0 }
        ')

    VCF_OTHER=$(bcftools view -H "$VCF" \
        | awk -F'\t' '
            {
                ref=$4
                alt=$5

                if (length(ref) == length(alt) && length(ref) == 1)
                    next

                # SNPs and simple indels have already been counted.
                # This category is deliberately conservative.
            }
        ')

else
    echo "WARNING: bcftools not found; VCF statistics will be unavailable."
fi

###############################################################################
# dnadiff
###############################################################################

DNADIFF_STATUS="NOT RUN"
DNADIFF_FILE="$TMPDIR_REPORT/dnadiff"

if command -v dnadiff >/dev/null 2>&1; then

    echo "Running dnadiff..."

    dnadiff \
        -p "$DNADIFF_FILE" \
        "$DRAFT" \
        "$POLISHED_FASTA" \
        >/dev/null 2>&1 || true

    DNADIFF_STATUS="COMPLETED"

else
    echo "WARNING: dnadiff not found; assembly comparison unavailable."
fi

###############################################################################
# QUAST
###############################################################################

QUAST_STATUS="NOT RUN"

if [[ -n "$QUAST_DIR" ]]; then

    if command -v quast.py >/dev/null 2>&1; then

        echo "Running QUAST..."

        mkdir -p "$QUAST_DIR"

        quast.py \
            "$DRAFT" \
            "$POLISHED_FASTA" \
            -o "$QUAST_DIR" \
            >/dev/null 2>&1 || true

        QUAST_STATUS="COMPLETED"

    else
        echo "WARNING: quast.py not found; QUAST not run."
        QUAST_STATUS="NOT AVAILABLE"
    fi
fi

###############################################################################
# Generate report
###############################################################################

{
    echo "============================================================"
    echo "DORADO POLISHING REPORT"
    echo "============================================================"
    echo
    echo "Generated: $(date)"
    echo
    echo "INPUT FILES"
    echo "------------------------------------------------------------"
    echo "Draft assembly:"
    echo "  $DRAFT"
    echo
    echo "Polished assembly:"
    echo "  $POLISHED"
    echo
    echo "Dorado VCF:"
    echo "  $VCF"
    echo
    echo
    echo "ASSEMBLY COMPARISON"
    echo "------------------------------------------------------------"
    printf "%-25s %15s %15s %15s\n" \
        "Metric" "Draft" "Polished" "Difference"
    echo "------------------------------------------------------------"

    printf "%-25s %15s %15s %15s\n" \
        "Assembly length (bp)" \
        "$DRAFT_LENGTH" \
        "$POLISHED_LENGTH" \
        "$LENGTH_DIFF"

    printf "%-25s %15s %15s %15s\n" \
        "Number of contigs" \
        "$DRAFT_CONTIGS" \
        "$POLISHED_CONTIGS" \
        "$((POLISHED_CONTIGS - DRAFT_CONTIGS))"

    printf "%-25s %15s %15s %15s\n" \
        "GC (%)" \
        "$DRAFT_GC" \
        "$POLISHED_GC" \
        "NA"

    echo
    echo
    echo "DORADO VARIANTS"
    echo "------------------------------------------------------------"
    echo "Total VCF records:       $VCF_TOTAL"
    echo "SNPs:                    $VCF_SNP"
    echo "Insertions:              $VCF_INS"
    echo "Deletions:               $VCF_DEL"
    echo
    echo
    echo "INTERPRETATION"
    echo "------------------------------------------------------------"
    echo "The VCF reports variants identified during Dorado polishing."
    echo "The assembly comparison reports sequence differences between"
    echo "the original draft and the polished consensus."
    echo
    echo "A difference between the assemblies does not by itself"
    echo "demonstrate that the polished assembly is more accurate."
    echo "Accuracy should ideally be assessed using an independent"
    echo "reference or independent validation data."
    echo
    echo
    echo "DNADIFF"
    echo "------------------------------------------------------------"
    echo "Status: $DNADIFF_STATUS"

    if [[ "$DNADIFF_STATUS" == "COMPLETED" && -f "${DNADIFF_FILE}.report" ]]; then
        echo
        echo "Key dnadiff results:"
        grep -E '^[[:space:]]*(NUCmer|1-to-1|AvgIdentity|SNPs|Indels)' \
            "${DNADIFF_FILE}.report" \
            || true
    fi

    echo
    echo
    echo "QUAST"
    echo "------------------------------------------------------------"
    echo "Status: $QUAST_STATUS"

    if [[ "$QUAST_STATUS" == "COMPLETED" ]]; then
        echo "Results directory:"
        echo "  $QUAST_DIR"
    fi

    echo
    echo "============================================================"
    echo "END OF REPORT"
    echo "============================================================"

} > "$REPORT"

###############################################################################
# Final message
###############################################################################

echo
echo "Report written to:"
echo "  $REPORT"
echo

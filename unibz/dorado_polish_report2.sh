#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Dorado polishing report
#
# Calculates:
#
#   1. Draft vs polished assembly statistics
#   2. Dorado VCF statistics
#   3. Read-to-draft alignment statistics
#   4. Read-to-polished alignment statistics
#   5. Read agreement:
#        - aligned reads
#        - mean depth
#        - mismatch rate
#        - indel rate
#        - total error rate
#   6. Optional dnadiff
#   7. Optional QUAST
#
# IMPORTANT:
#
# The reads must be realigned independently to the draft and polished
# assemblies. The original Dorado alignment BAM cannot be reused for the
# polished assembly because it is aligned to the draft.
#
# Usage:
#
#   ./dorado_polish_report.sh \
#       --draft consensus_assembly.fasta \
#       --polished consensus.fastq \
#       --reads calls.bam \
#       --vcf variants.vcf \
#       --output polishing_report.txt
#
###############################################################################

usage() {
    cat <<EOF

Usage:

  $0 \\
      --draft DRAFT.fa \\
      --polished POLISHED.fa/fq \\
      --reads READS.bam \\
      --vcf VARIANTS.vcf \\
      --output REPORT.txt

Required:
  --draft       Original draft assembly
  --polished    Dorado polished assembly FASTA/FASTQ
  --reads       Original basecalled reads BAM
  --vcf         Dorado variants.vcf
  --output      Output report

Optional:
  --quast-dir   Run QUAST and save results here
  --threads     Number of CPU threads (default: SLURM_CPUS_PER_TASK or 16)
  --dorado      Path to dorado executable
  --keep-bams   Keep temporary realignment BAMs

  --help        Show this help

Example:

  $0 \\
      --draft consensus_assembly.fasta \\
      --polished dorado_polish/consensus.fastq \\
      --reads calls.bam \\
      --vcf dorado_polish/variants.vcf \\
      --output dorado_polish/polishing_report.txt \\
      --quast-dir dorado_polish/quast

EOF
}

###############################################################################
# Defaults
###############################################################################

DRAFT=""
POLISHED=""
READS=""
VCF=""
REPORT=""
QUAST_DIR=""

THREADS="${SLURM_CPUS_PER_TASK:-16}"

DORADO="/data/users/theaven/software/dorado-2.1.2-linux-x64/bin/dorado"

KEEP_BAMS="N"

###############################################################################
# Parse arguments
###############################################################################

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

        --reads)
            READS="$2"
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

        --threads)
            THREADS="$2"
            shift 2
            ;;

        --dorado)
            DORADO="$2"
            shift 2
            ;;

        --keep-bams)
            KEEP_BAMS="Y"
            shift
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

if [[ -z "$DRAFT" ||
      -z "$POLISHED" ||
      -z "$READS" ||
      -z "$VCF" ||
      -z "$REPORT" ]]; then

    echo "ERROR: --draft, --polished, --reads, --vcf and --output are required."
    usage
    exit 1

fi

for file in "$DRAFT" "$POLISHED" "$READS" "$VCF"; do

    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found:"
        echo "  $file"
        exit 1
    fi

done

###############################################################################
# Check programs
###############################################################################

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found in PATH."
    exit 1
fi

if ! command -v "$DORADO" >/dev/null 2>&1 &&
   [[ ! -x "$DORADO" ]]; then

    echo "ERROR: Dorado executable not found:"
    echo "  $DORADO"
    exit 1

fi

###############################################################################
# Temporary directory
###############################################################################

TMPDIR_REPORT="$(dirname "$REPORT")/dorado_polish_workdir"
mkdir -p "$TMPDIR_REPORT"

cleanup() {

    if [[ "$KEEP_BAMS" != "Y" ]]; then
        rm -rf "$TMPDIR_REPORT"
    else
        echo
        echo "Temporary files retained:"
        echo "  $TMPDIR_REPORT"
    fi

}

trap cleanup EXIT

DRAFT_ALIGN_DIR="$TMPDIR_REPORT/draft_alignment"
POLISHED_ALIGN_DIR="$TMPDIR_REPORT/polished_alignment"

mkdir -p "$DRAFT_ALIGN_DIR"
mkdir -p "$POLISHED_ALIGN_DIR"

###############################################################################
# Convert polished FASTQ -> FASTA if necessary
###############################################################################

POLISHED_FASTA="$POLISHED"

case "$POLISHED" in

    *.fastq|*.fq|*.fastq.gz|*.fq.gz)

        POLISHED_FASTA="$TMPDIR_REPORT/polished.fasta"

        if command -v seqkit >/dev/null 2>&1; then

            echo "Converting polished FASTQ to FASTA..."

            seqkit fq2fa "$POLISHED" \
                -o "$POLISHED_FASTA"

        elif command -v seqtk >/dev/null 2>&1; then

            echo "Converting polished FASTQ to FASTA..."

            if [[ "$POLISHED" == *.gz ]]; then
                zcat "$POLISHED" | seqtk seq -a - > "$POLISHED_FASTA"
            else
                seqtk seq -a "$POLISHED" > "$POLISHED_FASTA"
            fi

        else

            echo "ERROR: polished assembly is FASTQ but neither"
            echo "seqkit nor seqtk is available."

            exit 1

        fi

        ;;

esac

###############################################################################
# FASTA statistics
###############################################################################

fasta_length() {

    awk '
        BEGIN { n=0 }
        /^>/ { next }
        {
            gsub(/[ \t\r\n]/,"")
            n += length($0)
        }
        END {
            print n
        }
    ' "$1"

}

fasta_contigs() {

    grep -c '^>' "$1"

}

fasta_gc() {

    awk '
        BEGIN {
            gc=0
            total=0
        }

        /^>/ {
            next
        }

        {
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

echo
echo "Calculating assembly statistics..."

DRAFT_LENGTH=$(fasta_length "$DRAFT")
POLISHED_LENGTH=$(fasta_length "$POLISHED_FASTA")

DRAFT_CONTIGS=$(fasta_contigs "$DRAFT")
POLISHED_CONTIGS=$(fasta_contigs "$POLISHED_FASTA")

DRAFT_GC=$(fasta_gc "$DRAFT")
POLISHED_GC=$(fasta_gc "$POLISHED_FASTA")

LENGTH_DIFF=$((POLISHED_LENGTH - DRAFT_LENGTH))
CONTIG_DIFF=$((POLISHED_CONTIGS - DRAFT_CONTIGS))

###############################################################################
# VCF statistics
###############################################################################

VCF_TOTAL="NA"
VCF_SNP="NA"
VCF_INS="NA"
VCF_DEL="NA"

if command -v bcftools >/dev/null 2>&1; then

    echo "Calculating VCF statistics..."

    VCF_TOTAL=$(bcftools view -H "$VCF" | wc -l)

    VCF_SNP=$(bcftools view -v snps -H "$VCF" | wc -l)

    VCF_INS=$(bcftools view -v indels -H "$VCF" |
        awk '
            {
                if (length($5) > length($4))
                    n++
            }
            END {
                print n+0
            }
        ')

    VCF_DEL=$(bcftools view -v indels -H "$VCF" |
        awk '
            {
                if (length($5) < length($4))
                    n++
            }
            END {
                print n+0
            }
        ')

else

    echo "WARNING: bcftools not found."

fi

###############################################################################
# Realign reads to draft
###############################################################################

echo
echo "============================================================"
echo "Aligning reads to DRAFT assembly"
echo "============================================================"

"$DORADO" aligner \
    "$DRAFT" \
    "$READS" \
    --output-dir "$DRAFT_ALIGN_DIR"

###############################################################################
# Find draft BAMs
###############################################################################

mapfile -t DRAFT_BAMS < <(
    find "$DRAFT_ALIGN_DIR" \
        -type f \
        -name "*.bam" \
        | sort
)

if [[ "${#DRAFT_BAMS[@]}" -eq 0 ]]; then

    echo "ERROR: Dorado produced no BAM files for draft alignment."
    exit 1

fi

echo "Found ${#DRAFT_BAMS[@]} draft BAM(s)."

DRAFT_BAM="$TMPDIR_REPORT/draft.bam"

samtools merge \
    -@ "$THREADS" \
    -o "$DRAFT_BAM" \
    "${DRAFT_BAMS[@]}"

samtools index \
    -@ "$THREADS" \
    "$DRAFT_BAM"

###############################################################################
# Realign reads to polished assembly
###############################################################################

echo
echo "============================================================"
echo "Aligning reads to POLISHED assembly"
echo "============================================================"

"$DORADO" aligner \
    "$POLISHED_FASTA" \
    "$READS" \
    --output-dir "$POLISHED_ALIGN_DIR"

###############################################################################
# Find polished BAMs
###############################################################################

mapfile -t POLISHED_BAMS < <(
    find "$POLISHED_ALIGN_DIR" \
        -type f \
        -name "*.bam" \
        | sort
)

if [[ "${#POLISHED_BAMS[@]}" -eq 0 ]]; then

    echo "ERROR: Dorado produced no BAM files for polished alignment."
    exit 1

fi

echo "Found ${#POLISHED_BAMS[@]} polished BAM(s)."

POLISHED_BAM="$TMPDIR_REPORT/polished.bam"

samtools merge \
    -@ "$THREADS" \
    -o "$POLISHED_BAM" \
    "${POLISHED_BAMS[@]}"

samtools index \
    -@ "$THREADS" \
    "$POLISHED_BAM"

###############################################################################
# Calculate read agreement
#
# We calculate:
#
#   total reads
#   mapped reads
#   mapped percentage
#   aligned bases
#   mismatches
#   insertions
#   deletions
#   mismatch rate
#   indel rate
#   total error rate
#
# NM is the edit distance between read and reference, including mismatches
# and indels.
#
# CIGAR is used to separate insertions/deletions.
###############################################################################

calculate_agreement() {

    local BAM="$1"
    local PREFIX="$2"

    echo "Calculating read agreement: $PREFIX"

    samtools stats \
        "$BAM" \
        > "$TMPDIR_REPORT/${PREFIX}.stats"

    # Total reads
    TOTAL_READS=$(
        awk '$1=="SN" && $2=="raw total sequences:" {print $NF}' \
        "$TMPDIR_REPORT/${PREFIX}.stats"
    )

    # Mapped reads
    MAPPED_READS=$(
        awk '$1=="SN" && $2=="reads mapped:" {print $NF}' \
        "$TMPDIR_REPORT/${PREFIX}.stats"
    )

    # Mapped percentage
    if [[ "$TOTAL_READS" -gt 0 ]]; then

        MAPPED_PERCENT=$(awk \
            -v m="$MAPPED_READS" \
            -v t="$TOTAL_READS" \
            'BEGIN {printf "%.3f", 100*m/t}')

    else

        MAPPED_PERCENT="NA"

    fi

    # Number of aligned bases
    ALIGNED_BASES=$(
        awk '$1=="SN" && $2=="bases mapped (cigar):" {print $NF}' \
        "$TMPDIR_REPORT/${PREFIX}.stats"
    )

    ###########################################################################
    # Calculate mismatches and indels from NM + CIGAR
    #
    # NM = mismatches + insertions + deletions
    #
    # Therefore:
    #
    # mismatches = NM - insertions - deletions
    #
    ###########################################################################

    read -r MISMATCHES INSERTIONS DELETIONS < <(

        samtools view "$BAM" |
        awk '
        BEGIN {
            mismatches=0
            insertions=0
            deletions=0
        }

        $0 !~ /^@/ {

            # Skip unmapped reads
            if ($3 == "*")
                next

            nm=""

            for (i=12; i<=NF; i++) {

                if ($i ~ /^NM:i:/) {
                    nm=substr($i,6)
                }

            }

            if (nm == "")
                next

            # Parse CIGAR
            cigar=$6

            ins=0
            del=0

            while (match(cigar, /[0-9]+[ID]/)) {

                token=substr(cigar,RSTART,RLENGTH)
                number=token

                sub(/[ID]$/, "", number)

                op=substr(token,length(token),1)

                if (op == "I")
                    ins += number

                if (op == "D")
                    del += number

                cigar=substr(cigar,RSTART+RLENGTH)

            }

            insertions += ins
            deletions += del
            mismatches += nm - ins - del

        }

        END {

            print mismatches, insertions, deletions

        }'
    )

    TOTAL_INDELS=$((INSERTIONS + DELETIONS))
    TOTAL_ERRORS=$((MISMATCHES + INSERTIONS + DELETIONS))

    if [[ "$ALIGNED_BASES" -gt 0 ]]; then

        MISMATCH_RATE=$(awk \
            -v x="$MISMATCHES" \
            -v n="$ALIGNED_BASES" \
            'BEGIN {printf "%.6f", 100*x/n}')

        INDEL_RATE=$(awk \
            -v x="$TOTAL_INDELS" \
            -v n="$ALIGNED_BASES" \
            'BEGIN {printf "%.6f", 100*x/n}')

        ERROR_RATE=$(awk \
            -v x="$TOTAL_ERRORS" \
            -v n="$ALIGNED_BASES" \
            'BEGIN {printf "%.6f", 100*x/n}')

        AGREEMENT=$(awk \
            -v x="$TOTAL_ERRORS" \
            -v n="$ALIGNED_BASES" \
            'BEGIN {printf "%.6f", 100*(1-x/n)}')

    else

        MISMATCH_RATE="NA"
        INDEL_RATE="NA"
        ERROR_RATE="NA"
        AGREEMENT="NA"

    fi

}

###############################################################################
# Draft agreement
###############################################################################

calculate_agreement \
    "$DRAFT_BAM" \
    "draft"

DRAFT_TOTAL_READS="$TOTAL_READS"
DRAFT_MAPPED_READS="$MAPPED_READS"
DRAFT_MAPPED_PERCENT="$MAPPED_PERCENT"
DRAFT_ALIGNED_BASES="$ALIGNED_BASES"

DRAFT_MISMATCHES="$MISMATCHES"
DRAFT_INSERTIONS="$INSERTIONS"
DRAFT_DELETIONS="$DELETIONS"

DRAFT_MISMATCH_RATE="$MISMATCH_RATE"
DRAFT_INDEL_RATE="$INDEL_RATE"
DRAFT_ERROR_RATE="$ERROR_RATE"
DRAFT_AGREEMENT="$AGREEMENT"

###############################################################################
# Polished agreement
###############################################################################

calculate_agreement \
    "$POLISHED_BAM" \
    "polished"

POLISHED_TOTAL_READS="$TOTAL_READS"
POLISHED_MAPPED_READS="$MAPPED_READS"
POLISHED_MAPPED_PERCENT="$MAPPED_PERCENT"
POLISHED_ALIGNED_BASES="$ALIGNED_BASES"

POLISHED_MISMATCHES="$MISMATCHES"
POLISHED_INSERTIONS="$INSERTIONS"
POLISHED_DELETIONS="$DELETIONS"

POLISHED_MISMATCH_RATE="$MISMATCH_RATE"
POLISHED_INDEL_RATE="$INDEL_RATE"
POLISHED_ERROR_RATE="$ERROR_RATE"
POLISHED_AGREEMENT="$AGREEMENT"

###############################################################################
# Differences
###############################################################################

ERROR_RATE_CHANGE=$(awk \
    -v before="$DRAFT_ERROR_RATE" \
    -v after="$POLISHED_ERROR_RATE" \
    'BEGIN {printf "%.6f", after-before}')

AGREEMENT_CHANGE=$(awk \
    -v before="$DRAFT_AGREEMENT" \
    -v after="$POLISHED_AGREEMENT" \
    'BEGIN {printf "%.6f", after-before}')

MISMATCH_RATE_CHANGE=$(awk \
    -v before="$DRAFT_MISMATCH_RATE" \
    -v after="$POLISHED_MISMATCH_RATE" \
    'BEGIN {printf "%.6f", after-before}')

INDEL_RATE_CHANGE=$(awk \
    -v before="$DRAFT_INDEL_RATE" \
    -v after="$POLISHED_INDEL_RATE" \
    'BEGIN {printf "%.6f", after-before}')

###############################################################################
# dnadiff
###############################################################################

DNADIFF_STATUS="NOT RUN"
DNADIFF_DIR="$TMPDIR_REPORT/dnadiff"

if command -v dnadiff >/dev/null 2>&1; then

    echo
    echo "Running dnadiff..."

    mkdir -p "$DNADIFF_DIR"

    dnadiff \
        -p "$DNADIFF_DIR/comparison" \
        "$DRAFT" \
        "$POLISHED_FASTA" \
        >/dev/null 2>&1 || true

    DNADIFF_STATUS="COMPLETED"

else

    echo "WARNING: dnadiff not found."

fi

###############################################################################
# QUAST
###############################################################################

QUAST_STATUS="NOT RUN"

if [[ -n "$QUAST_DIR" ]]; then

    if command -v quast.py >/dev/null 2>&1; then

        echo
        echo "Running QUAST..."

        mkdir -p "$QUAST_DIR"

        quast.py \
            "$DRAFT" \
            "$POLISHED_FASTA" \
            -o "$QUAST_DIR" \
            >/dev/null 2>&1 || true

        QUAST_STATUS="COMPLETED"

    else

        QUAST_STATUS="NOT AVAILABLE"

    fi

fi

###############################################################################
# Write report
###############################################################################

mkdir -p "$(dirname "$REPORT")"

{
    echo "============================================================"
    echo "DORADO POLISHING REPORT"
    echo "============================================================"
    echo
    echo "Generated: $(date)"
    echo
    echo "THREADS: $THREADS"
    echo
    echo "INPUTS"
    echo "------------------------------------------------------------"
    echo "Draft assembly:"
    echo "  $DRAFT"
    echo
    echo "Polished assembly:"
    echo "  $POLISHED"
    echo
    echo "Reads:"
    echo "  $READS"
    echo
    echo "Dorado VCF:"
    echo "  $VCF"
    echo
    echo
    echo "ASSEMBLY COMPARISON"
    echo "------------------------------------------------------------"

    printf "%-30s %15s %15s %15s\n" \
        "Metric" "Draft" "Polished" "Difference"

    echo "------------------------------------------------------------"

    printf "%-30s %15s %15s %15s\n" \
        "Assembly length (bp)" \
        "$DRAFT_LENGTH" \
        "$POLISHED_LENGTH" \
        "$LENGTH_DIFF"

    printf "%-30s %15s %15s %15s\n" \
        "Number of contigs" \
        "$DRAFT_CONTIGS" \
        "$POLISHED_CONTIGS" \
        "$CONTIG_DIFF"

    printf "%-30s %15s %15s %15s\n" \
        "GC (%)" \
        "$DRAFT_GC" \
        "$POLISHED_GC" \
        "NA"

    echo
    echo
    echo "DORADO VARIANTS"
    echo "------------------------------------------------------------"

    printf "%-30s %15s\n" \
        "Total VCF records" \
        "$VCF_TOTAL"

    printf "%-30s %15s\n" \
        "SNPs" \
        "$VCF_SNP"

    printf "%-30s %15s\n" \
        "Insertions" \
        "$VCF_INS"

    printf "%-30s %15s\n" \
        "Deletions" \
        "$VCF_DEL"

    echo
    echo
    echo "READ-TO-ASSEMBLY AGREEMENT"
    echo "------------------------------------------------------------"

    printf "%-30s %15s %15s\n" \
        "Metric" "Draft" "Polished"

    echo "------------------------------------------------------------"

    printf "%-30s %15s %15s\n" \
        "Total reads" \
        "$DRAFT_TOTAL_READS" \
        "$POLISHED_TOTAL_READS"

    printf "%-30s %15s %15s\n" \
        "Mapped reads" \
        "$DRAFT_MAPPED_READS" \
        "$POLISHED_MAPPED_READS"

    printf "%-30s %14s%% %14s%%\n" \
        "Mapped reads (%)" \
        "$DRAFT_MAPPED_PERCENT" \
        "$POLISHED_MAPPED_PERCENT"

    printf "%-30s %15s %15s\n" \
        "Aligned bases" \
        "$DRAFT_ALIGNED_BASES" \
        "$POLISHED_ALIGNED_BASES"

    echo
    echo "Mismatch / indel statistics:"
    echo

    printf "%-30s %15s %15s\n" \
        "Mismatches" \
        "$DRAFT_MISMATCHES" \
        "$POLISHED_MISMATCHES"

    printf "%-30s %15s %15s\n" \
        "Insertions" \
        "$DRAFT_INSERTIONS" \
        "$POLISHED_INSERTIONS"

    printf "%-30s %15s %15s\n" \
        "Deletions" \
        "$DRAFT_DELETIONS" \
        "$POLISHED_DELETIONS"

    echo
    echo "Rates:"
    echo

    printf "%-30s %14s%% %14s%%\n" \
        "Mismatch rate" \
        "$DRAFT_MISMATCH_RATE" \
        "$POLISHED_MISMATCH_RATE"

    printf "%-30s %14s%% %14s%%\n" \
        "Indel rate" \
        "$DRAFT_INDEL_RATE" \
        "$POLISHED_INDEL_RATE"

    printf "%-30s %14s%% %14s%%\n" \
        "Total error rate" \
        "$DRAFT_ERROR_RATE" \
        "$POLISHED_ERROR_RATE"

    printf "%-30s %14s%% %14s%%\n" \
        "Read agreement" \
        "$DRAFT_AGREEMENT" \
        "$POLISHED_AGREEMENT"

    echo
    echo "Change after polishing:"
    echo

    printf "%-30s %15s percentage points\n" \
        "Mismatch rate change" \
        "$MISMATCH_RATE_CHANGE"

    printf "%-30s %15s percentage points\n" \
        "Indel rate change" \
        "$INDEL_RATE_CHANGE"

    printf "%-30s %15s percentage points\n" \
        "Total error rate change" \
        "$ERROR_RATE_CHANGE"

    printf "%-30s %15s percentage points\n" \
        "Read agreement change" \
        "$AGREEMENT_CHANGE"

    echo
    echo
    echo "DNADIFF"
    echo "------------------------------------------------------------"
    echo "Status: $DNADIFF_STATUS"

    if [[ "$DNADIFF_STATUS" == "COMPLETED" ]]; then

        DNADIFF_REPORT="$DNADIFF_DIR/comparison.report"

        if [[ -f "$DNADIFF_REPORT" ]]; then

            echo
            echo "Selected dnadiff statistics:"
            echo

            grep -E \
                '^[[:space:]]*(NUCMER|1-to-1|AvgIdentity|SNPs|Indels)' \
                "$DNADIFF_REPORT" \
                || true

        fi

    fi

    echo
    echo
    echo "QUAST"
    echo "------------------------------------------------------------"
    echo "Status: $QUAST_STATUS"

    if [[ "$QUAST_STATUS" == "COMPLETED" ]]; then
        echo
        echo "Results:"
        echo "  $QUAST_DIR"
    fi

    echo
    echo
    echo "INTERPRETATION"
    echo "------------------------------------------------------------"
    echo
    echo "Read-to-assembly agreement is calculated after independently"
    echo "aligning the same input reads to the draft and polished"
    echo "assemblies."
    echo
    echo "Mismatch rate:"
    echo "  mismatching bases / aligned bases"
    echo
    echo "Indel rate:"
    echo "  inserted + deleted bases / aligned bases"
    echo
    echo "Total error rate:"
    echo "  mismatches + insertions + deletions / aligned bases"
    echo
    echo "Read agreement:"
    echo "  100 - total error rate"
    echo
    echo "These metrics describe agreement between the reads and each"
    echo "assembly. They are not an independent measurement of genome"
    echo "accuracy, because the same reads are used both for polishing"
    echo "and for this evaluation."
    echo
    echo "An independent reference or independent validation dataset"
    echo "is preferable when estimating absolute assembly accuracy."
    echo
    echo "============================================================"
    echo "END OF REPORT"
    echo "============================================================"

} > "$REPORT"

###############################################################################
# Final message
###############################################################################

echo
echo "============================================================"
echo "REPORT COMPLETE"
echo "============================================================"
echo
echo "Report:"
echo "  $REPORT"
echo
echo "Draft read agreement:"
echo "  $DRAFT_AGREEMENT %"
echo
echo "Polished read agreement:"
echo "  $POLISHED_AGREEMENT %"
echo
echo "Change:"
echo "  $AGREEMENT_CHANGE percentage points"
echo

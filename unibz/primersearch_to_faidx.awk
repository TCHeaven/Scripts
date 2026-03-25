# Usage:
#   samtools faidx genome.fasta
#   awk -v GENOME="genome.fasta" -f primersearch_to_faidx.awk genome.fasta.fai primer_results.txt > extract_cmds.sh
#   bash extract_cmds.sh > amplicons.fasta

BEGIN { FS = "[\t ]+" }

# Read contig lengths from the .fai (1st file)
NR==FNR {
  seqlen[$1] = $2
  next
}

# Primer name ITS1-ITS4_Pyt
/^Primer[ \t]+name[ \t]+/ {
  primer = $3
  next
}

# Amplimer 1
/^Amplimer[ \t]+[0-9]+/ {
  amp = $2
  contig = ""
  fpos = 0
  alen = 0
  next
}

#         Sequence: JADDUH010000185.1
/^[ \t]*Sequence:[ \t]*/ {
  # capture the first non-whitespace token after "Sequence:"
  if (match($0, /^[ \t]*Sequence:[ \t]*([^ \t\r\n]+)/, m)) {
    contig = m[1]
    gsub(/:$/, "", contig)   # safety: remove trailing colon if any
  }
  next
}

#         ... hits forward strand at 32872 with 1 mismatches
/hits forward strand at/ {
  if (match($0, / at ([0-9]+)/, m)) fpos = m[1] + 0
  next
}

#         Amplimer length: 2045 bp
/Amplimer length:/ {
  if (match($0, /Amplimer length:[ \t]*([0-9]+)/, m)) alen = m[1] + 0

  if (primer == "" || contig == "" || fpos == 0 || alen == 0) next

  # Optional sanity check: contig must exist in FASTA index
  if (!(contig in seqlen)) {
    # print a warning to stderr and skip
    printf("WARN: contig '%s' not found in .fai (primer=%s amplimer=%s)\n", contig, primer, amp) > "/dev/stderr"
    next
  }

  start = fpos
  end   = fpos + alen - 1
  if (end > seqlen[contig]) end = seqlen[contig]  # clamp just in case

  outname = primer "_Amplimer" amp "_" contig ":" start "-" end

  printf("samtools faidx %s %s:%d-%d | sed '1s/^>.*/>%s/'\n",
         GENOME, contig, start, end, outname)
  next
}

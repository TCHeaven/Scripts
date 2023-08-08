import sys
import os
import pysam
from pyfaidx import Fasta

multi_fasta_file = sys.argv[1]
gff3_file = sys.argv[2]
output_file = sys.argv[3]

def splice_genes(gff3_file, multi_fasta_file, output_file):
    # Parse the GFF3 file to extract gene coordinates
    gene_coordinates = {}
    CDS_coordinates = {}
    with open(gff3_file, "r") as gff_file:
        #print("gff open")
        i = 1
        for line in gff_file:
            #print("processing gff")
            if line.startswith("#"):
                continue
            gene_info = line.strip().split("\t")
            if gene_info[2] == "gene":
                gene_id = gene_info[8].split(";")[0].replace("ID=", "")
                gene_chrom = gene_info[0]  # Store the chromosome name from the GFF3 file
                gene_start = int(gene_info[3]) - 1
                print("Poaition before gene start:")
                print(gene_start)
                gene_end = int(gene_info[4])
            if gene_info[2] == "CDS":
                CDS_id = gene_info[8].split(";")[0].replace("ID=", "") + str(i)
                i += 1
                CDS_chrom = gene_info[0]  # Store the chromosome name from the GFF3 file
                CDS_start = int(gene_info[3]) - gene_start
                print("CDS start position relative to gene:")
                print(CDS_start)
                CDS_end = int(gene_info[4]) - gene_start
                CDS_coordinates[CDS_id] = (CDS_chrom, CDS_start, CDS_end)
    with open(multi_fasta_file, "r") as fasta_file, open(output_file, "w") as output:
        for fastaline in fasta_file:
            fastaline = fastaline.strip()
            if fastaline.startswith(">"):
                current_seq_id = fastaline[1:]
                output.write(f"\n>{current_seq_id}_CDS\n")
                print(f"\n{current_seq_id}")
                current_seq = ""
            else:
                #print("splicing")
                current_seq += fastaline  # Script will only work if fasta is in single line format 
                #print(current_seq)
                for cds_id, (CDS_chrom, CDS_start, CDS_end) in sorted(CDS_coordinates.items()):
                    if CDS_start <= CDS_end <= len(current_seq):
                        print(cds_id)
                        start_pos = CDS_start - 1
                        end_pos = CDS_end
                        print(current_seq[start_pos:end_pos])
                        output.write(current_seq[start_pos:end_pos])
                    else:
                        print("ERROR")
                        print(cds_id)
                        print("CDS start:")
                        print(CDS_start)
                        print("gene length:")
                        print(len(current_seq))
                        print("CDS end:")
                        print(CDS_end)

splice_genes(gff3_file, multi_fasta_file, output_file)
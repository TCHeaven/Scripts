import sys
import os 
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.AlignInfo import SummaryInfo

def calculate_dnds_ratio(dn, ds):
    if dn != 0:
        return ds / dn
    else:
        return float("inf")

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py multifasta.fasta")
        sys.exit(1)
    multifasta = sys.argv[1]
    # Perform multiple sequence alignment
    clustalw_cline = ClustalwCommandline("/nbi/software/testing/bin/clustalw-2.1", infile=multifasta)
    clustalw_cline()
    alignment_file = os.path.abspath("temp.aln")
    alignment = AlignIO.read(alignment_file, "clustal")
    # Calculate dN and dS values for the whole alignment
    summary = SummaryInfo(alignment)
    dn = summary.difference_matrix(["A", "C", "G", "T"], output="ignore_missing")
    ds = summary.codon_usage_table().synonymous_matrix(["A", "C", "G", "T"])
    total_dn = sum(dn.values())
    total_ds = sum(ds.values())
    dnds_ratio = calculate_dnds_ratio(total_dn, total_ds)
    print(f"Total dN: {total_dn}, Total dS: {total_ds}, Overall dN/dS ratio: {dnds_ratio}")

if __name__ == "__main__":
    main()

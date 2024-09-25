import argparse
import pysam

def extract_af(vcf_file):
    af_dict = {}
    vcf = pysam.VariantFile(vcf_file, "r")
    for rec in vcf.fetch():
        chrom = rec.chrom
        pos = rec.pos
        info = rec.info
        af = info.get("AF", None)
        if af is not None:
            if isinstance(af, tuple):
                af = af[0]
            elif isinstance(af, list) and len(af) > 0:
                af = af[0] 
            af_dict[(chrom, pos)] = af
    return af_dict

def main(group1_vcf, group2_vcf, output_file):
    af_group1 = extract_af(group1_vcf)
    af_group2 = extract_af(group2_vcf)
    with open(output_file, "w") as out:
        out.write("CHROM\tPOS\tGROUP1_AF\tGROUP2_AF\tDIFFERENCE\n")
        all_positions = set(af_group1.keys()).union(set(af_group2.keys()))
        for pos in all_positions:
            af1 = af_group1.get(pos, "NA")
            af2 = af_group2.get(pos, "NA")
            if af1 != "NA":
                af1 = float(af1)
            if af2 != "NA":
                af2 = float(af2)
            if af1 != "NA" and af2 != "NA":
                difference = abs(af2 - af1)
            else:
                difference = "NA"
            out.write(f"{pos[0]}\t{pos[1]}\t{af1}\t{af2}\t{difference}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare AF values between two VCF files.")
    parser.add_argument("group1_vcf", help="Path to the first group VCF file")
    parser.add_argument("group2_vcf", help="Path to the second group VCF file")
    parser.add_argument("output_file", help="Path to the output file")
    args = parser.parse_args()
    main(args.group1_vcf, args.group2_vcf, args.output_file)

import argparse
import pysam
import multiprocessing as mp

def extract_af_chunk(vcf_file, chunk_start, chunk_end):
    af_dict = {}
    vcf = pysam.VariantFile(vcf_file, "r")
    current_position = 0
    print(f"Processing chunk from {chunk_start} to {chunk_end}")
    for rec in vcf.fetch():
        if current_position >= chunk_end:
            break
        if current_position >= chunk_start:
            chrom = rec.chrom
            pos = rec.pos
            info = rec.info
            af = info.get("AF", None)
            if af is not None:
                if isinstance(af, tuple):
                    af = af[0]
                elif isinstance(af, list) and len(af) > 0:
                    af = af[0]
                print(f"Record: {chrom} {pos}, AF: {af}")  
                af_dict[(chrom, pos)] = af
        current_position += 1
    return af_dict

def extract_af(vcf_file, chunk_size, num_chunks):
    vcf = pysam.VariantFile(vcf_file, "r")
    num_records = sum(1 for _ in vcf.fetch())
    print(f"Total number of records: {num_records}")
    chunk_ranges = [(i * chunk_size, min((i + 1) * chunk_size, num_records)) for i in range(num_chunks)]
    with mp.Pool(processes=num_chunks) as pool:
        results = pool.starmap(extract_af_chunk, [(vcf_file, start, end) for start, end in chunk_ranges])
    af_dict = {}
    for result in results:
        af_dict.update(result)
    return af_dict

def main(group1_vcf, group2_vcf, output_file, chunk_size=10000, num_chunks=8):
    af_group1 = extract_af(group1_vcf, chunk_size, num_chunks)
    af_group2 = extract_af(group2_vcf, chunk_size, num_chunks)
    with open(output_file, "w") as out:
        out.write("CHROM\tPOS\tGROUP1_AF\tGROUP2_AF\tDIFFERENCE\n")
        buffer = []
        all_positions = set(af_group1.keys()).union(set(af_group2.keys()))
        for pos in all_positions:
            af1 = af_group1.get(pos, "NA")
            af2 = af_group2.get(pos, "NA")
            try:
                af1 = float(af1) if af1 != "NA" else "NA"
            except ValueError:
                af1 = "NA"
            try:
                af2 = float(af2) if af2 != "NA" else "NA"
            except ValueError:
                af2 = "NA"
            if af1 != "NA" and af2 != "NA":
                difference = abs(af2 - af1)
            else:
                difference = "NA"
            buffer.append(f"{pos[0]}\t{pos[1]}\t{af1}\t{af2}\t{difference}\n")
            if len(buffer) >= 10000:
                out.writelines(buffer)
                buffer = []
        if buffer:
            out.writelines(buffer)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare AF values between two VCF files.")
    parser.add_argument("group1_vcf", help="Path to the first group VCF file")
    parser.add_argument("group2_vcf", help="Path to the second group VCF file")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("--chunk_size", type=int, default=10000, help="Size of each chunk (default: 10000)")
    parser.add_argument("--num_chunks", type=int, default=8, help="Number of parallel chunks (default: 8)")
    args = parser.parse_args()
    main(args.group1_vcf, args.group2_vcf, args.output_file, args.chunk_size, args.num_chunks)

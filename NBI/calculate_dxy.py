import sys
import pysam
import numpy as np

def compute_dxy(vcf_file, pop1_samples, pop2_samples):
    try:
        vcf = pysam.VariantFile(vcf_file, 'r')
    except Exception as e:
        print(f"Error opening VCF file: {e}")
        sys.exit(1)
    dxy_values = []
    for rec in vcf.fetch():
        # Extract genotypes for the two populations
        pop1_genotypes = [rec.samples[sample]['GT'] for sample in pop1_samples if sample in rec.samples]
        pop2_genotypes = [rec.samples[sample]['GT'] for sample in pop2_samples if sample in rec.samples]
        if len(pop1_genotypes) == 0 or len(pop2_genotypes) == 0:
            continue
        # Compute pairwise differences
        pairwise_differences = []
        for gt1 in pop1_genotypes:
            for gt2 in pop2_genotypes:
                if None in gt1 or None in gt2:
                    continue
                diff = np.sum(np.array(gt1) != np.array(gt2))
                pairwise_differences.append(diff)
        if pairwise_differences:
            dxy_values.append(np.mean(pairwise_differences))
    return np.mean(dxy_values) if dxy_values else float('nan')

def read_sample_ids(file_path):
    try:
        with open(file_path, 'r') as file:
            sample_ids = [line.strip() for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading sample file {file_path}: {e}")
        sys.exit(1)
    return sample_ids

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <vcf_file> <pop1_samples> <pop2_samples>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    group1 = sys.argv[2]
    group2 = sys.argv[3]
    pop1_samples = read_sample_ids(group1)
    pop2_samples = read_sample_ids(group2)
    dxy = compute_dxy(vcf_file, pop1_samples, pop2_samples)
    print(f"Dxy: {dxy}")

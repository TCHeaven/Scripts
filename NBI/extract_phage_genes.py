import json
import sys

json_file = sys.argv[1]
fasta = sys.argv[2]
gff = sys.argv[3]
faa = sys.argv[4]
ID = sys.argv[5]
output_file = f"{ID}_phage_genes.txt"

def read_fasta_and_create_dict(fasta, ID):
    header_dict = {}
    id_counter = 1
    try:
        with open(fasta, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    header = line[1:].strip()  # Remove '>' and whitespace
                    first_item = header.split(' ')[0]
                    unique_id = f"gnl|JIC|{ID}_{id_counter}"
                    header_dict[unique_id] = first_item
                    id_counter += 1
    except FileNotFoundError:
        print(f"Error: The file '{fasta}' was not found.")
        sys.exit(1)
    except IOError as e:
        print(f"Error: An I/O error occurred: {e}")
        sys.exit(1)
    return header_dict

header_dict = read_fasta_and_create_dict(fasta, ID)

def load_gff(gff, header_dict):
    genes = []
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:  # Check if parts has enough elements
                continue 
            if parts[2] == "gene":
                start = int(parts[3])
                stop = int(parts[4])
                contig = header_dict.get(parts[0], "Unknown")
                attributes = parts[8]
                gene_id = attributes.split(';')[0].split('=')[1]
                gene_id = gene_id.replace('_gene', '')
                genes.append((start, stop, contig, gene_id))
    return genes

genes = load_gff(gff, header_dict)

# Load JSON data from file
def load_json(json_file):
    with open(json_file, 'r') as f:
        return json.load(f)

json_data = load_json(json_file)

def extract_genes_in_region(json_data, genes):
    extracted_genes = []
    for region in json_data:
        if "contig_start" in region:
            region_start = region["start"] - region["contig_start"] + 1
            region_stop = region["stop"] - region["contig_start"] + 1
            contig_info = region["contig_tag"]
            region_contig = contig_info.split(',')[0]
            for gene in genes:
                start, stop, contig, gene_id = gene
                if region_contig == contig:
                  if region_start <= start <= region_stop or region_start <= stop <= region_stop:
                      extracted_genes.append(gene_id)
        else:
            region_start = region["start"]
            region_stop = region["stop"]
            for gene in genes:
                start, stop, contig, gene_id = gene
                if region_start <= start <= region_stop or region_start <= stop <= region_stop:
                    extracted_genes.append(gene_id)
    return extracted_genes

phage_genes = extract_genes_in_region(json_data, genes)

def write_to_file(phage_genes, output_file):
    with open(output_file, 'w') as f:
        for gene_id in phage_genes:
            f.write(f"{gene_id}\n")

write_to_file(phage_genes, output_file)

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--id_file', required=True,
                        help='file path of removable ids')
    parser.add_argument('--input', required=True,
                        help='path of fasta sequence')
    parser.add_argument('--output', required=True,
                        help='path of selected fasta sequence')

    args = parser.parse_args()

    with open(args.id_file, mode='r', encoding='utf-8') as file:
        ids = set(file.read().splitlines())

    records = [SeqRecord(record.seq, id=record.id, description='', name='') for record in SeqIO.parse(args.input, 'fasta')
               if record.id not in ids]
    SeqIO.write(records, args.output, "fasta")

#!/usr/bin/env python
import pandas as pd
import argparse

def extract_flanking_sequences(table_file, fasta_file, output_file):
    # Read the table file into a pandas DataFrame
    df = pd.read_csv(table_file, sep='\t')

    # Read the FASTA file into a dictionary
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chrom = line.strip()[1:]
                chrom = chrom[0:8]+"_"+chrom[-8::]
                fasta_dict[chrom] = ''
            else:
                fasta_dict[chrom] += line.strip()

    # Extract the flanking sequences
    with open(output_file, 'w') as f:
        for index, row in df.iterrows():
            chrom = row['Chrom (or Genes)']
            start = int(row['Start'])
            stop = int(row['Stop'])

            # Extract the upstream and downstream sequences
            upstream_start = max(0, start - 100)
            downstream_stop = min(len(fasta_dict[chrom]), stop + 100)
            upstream_seq = fasta_dict[chrom][upstream_start:start].lower()
            downstream_seq = fasta_dict[chrom][stop:downstream_stop].lower()
            amplicon = fasta_dict[chrom][start:stop]
            outseq = f'{upstream_seq}{amplicon}{downstream_seq}'

            # Write the results to the output file
            f.write(f'>{chrom}!!AmpliconCoords={upstream_start}:{start}_{stop}:{downstream_stop}\n{outseq}\n\n')

def main():
    parser = argparse.ArgumentParser(description='Extract flanking sequences from a FASTA file based on a table of coordinates.')
    parser.add_argument('-t', '--table', required=True, help='Input table file containing chromosome coordinates.')
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file containing chromosome sequences.')
    parser.add_argument('-o', '--output', required=True, help='Output file for flanking sequences.')
    args = parser.parse_args()

    extract_flanking_sequences(args.table, args.fasta, args.output)

if __name__ == '__main__':
    main()

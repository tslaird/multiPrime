#!/usr/bin/env python
import csv
import argparse

def replace_headers(mapping_file, input_fasta, output_fasta):
    """
    Replaces headers in a FASTA file based on a mapping CSV file.

    Args:
        mapping_file (str): Path to the mapping CSV file.
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.
    """

    # Load the mapping from the CSV file
    mapping = {}
    with open(mapping_file, 'r') as f_map:
        reader = csv.reader(f_map)
        next(reader)  # Skip the header row
        for row in reader:
            OriginalHeader,HashHeader,CondensedHashHeader = row
            mapping[CondensedHashHeader] = OriginalHeader

    # Replace headers in the FASTA file
    with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                old_header = line.strip()[1:]
                old_header = old_header.split("!!")[0]
                new_header = mapping.get(old_header)
                if new_header:
                    f_out.write(f'>{new_header}\n')
                else:
                    print(f"Warning: No mapping found for header '{old_header}'")
                    f_out.write(line)
            else:
                f_out.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Replace headers in a FASTA file based on a mapping CSV file.')
    parser.add_argument('-m', '--mapping', required=True, help='Path to the mapping CSV file.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input FASTA file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output FASTA file.')
    args = parser.parse_args()

    replace_headers(args.mapping, args.input, args.output)
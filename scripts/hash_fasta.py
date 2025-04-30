#!/usr/bin/env python
import hashlib
import argparse
import csv

def hash_fasta(input_file, output_file, mapping_file, hash_algorithm='sha256'):
    """
    Maps FASTA sequences to a hash value and generates a new FASTA file with the hash value as the header.
    Also generates a table mapping the old headers to the new headers.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
        mapping_file (str): Path to the output mapping file.
        hash_algorithm (str, optional): Hash algorithm to use. Defaults to 'sha256'.

    Raises:
        ValueError: If the hash algorithm is not supported.
    """

    # Check if the hash algorithm is supported
    try:
        hashlib.new(hash_algorithm)
    except ValueError:
        raise ValueError(f"Unsupported hash algorithm: {hash_algorithm}")

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out, open(mapping_file, 'w', newline='') as f_map:
        writer = csv.writer(f_map)
        writer.writerow(['OriginalHeader', 'HashHeader','CondensedHashHeader'])  # Header row
        seq = ''
        header = ''
        for line in f_in:
            if line.startswith('>'):
                if seq:
                    # Hash the sequence and the old header
                    seq_hash = hashlib.new(hash_algorithm, seq.encode()).hexdigest()
                    header_hash = hashlib.new(hash_algorithm, header.encode()).hexdigest()
                    combined_hash = hashlib.new(hash_algorithm, (seq_hash + header_hash).encode()).hexdigest()
                    # Write the previous sequence to the output file
                    f_out.write(f'>{combined_hash}\n')
                    f_out.write(seq + '\n\n')
                    # Write the mapping to the mapping file
                    condensed_hash = combined_hash[0:8]+"_"+combined_hash[-8::]
                    writer.writerow([header[1:], combined_hash, condensed_hash])
                # Read the new header and reset the sequence
                header = line.strip()
                seq = ''
            else:
                seq += line.strip()
        # Hash the last sequence and the old header
        if seq:
            seq_hash = hashlib.new(hash_algorithm, seq.encode()).hexdigest()
            header_hash = hashlib.new(hash_algorithm, header.encode()).hexdigest()
            combined_hash = hashlib.new(hash_algorithm, (seq_hash + header_hash).encode()).hexdigest()
            # Write the last sequence to the output file
            f_out.write(f'>{combined_hash}\n')
            f_out.write(seq + '\n')
            # Write the mapping to the mapping file
            writer.writerow([header[1:], combined_hash, seq_hash])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hash FASTA sequences and generate a new FASTA file with the hash value as the header.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input FASTA file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output FASTA file.')
    parser.add_argument('-m', '--mapping', required=True, help='Path to the output mapping file.')
    parser.add_argument('-a', '--algorithm', default='sha256', help='Hash algorithm to use (default: sha256).')
    args = parser.parse_args()

    hash_fasta(args.input, args.output, args.mapping, args.algorithm)
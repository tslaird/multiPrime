#!/usr/bin/env python

import csv
import argparse

def load_fasta(fasta_file):
    seqs = {}
    with open(fasta_file, 'r') as f:
        seq = ''
        header = ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    seqs[header] = seq
                header = line.strip()[1:]
                seq = ''
            else:
                seq += line.strip()
        if seq:
            seqs[header] = seq
    return seqs

def match_nucleotides(n1, n2):
    # Degenerate nucleotide alphabet
    degenerate_alphabet = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }

    return n2.upper() in degenerate_alphabet.get(n1.upper(), [])


def match_primers_and_amplicons(table_file, primer_fasta, genome_fasta,hashmap_file, output_file):
    # Load the primer sequences from the primer FASTA file
    primer_seqs = load_fasta(primer_fasta)

    # Load the genome sequences from the genome FASTA file
    genome_seqs = load_fasta(genome_fasta)

    hashmap = {}
    hashmap = {}
    with open(hashmap_file, 'r') as f_map:
        reader = csv.reader(f_map)
        next(reader)  # Skip header row
        for row in reader:
            original_header, hash_header, condensed_hash_header = row
            hashmap[condensed_hash_header] = original_header
    # Read the table file and match the primer sequences and amplicon sequences
    with open(table_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
        reader = csv.reader(f_in, delimiter='\t')
        headers = next(reader)
        headers.extend(['Primer_F_seq', 'Primer_R_seq', 'AmpliconSequence', 'Primer_F_seq_mismasked','Primer_R_seq_mismasked', 'Primer_F_seq_mismatch', 'Primer_R_seq_mismatch','OriginalHeader'])
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(headers)
        for row in reader:
            chrom, start, stop, primer_f, primer_r, *rest = row
            primer_f_seq = primer_seqs.get(primer_f, '')
            primer_r_seq = primer_seqs.get(primer_r, '')
            amplicon_seq = ''
            if chrom in genome_seqs:
                start, stop = int(start), int(stop)
                amplicon_seq = genome_seqs[chrom][start:stop+len(primer_r_seq)]
            primer_f_seq_aligned = align_primer(primer_f_seq, amplicon_seq[:len(primer_f_seq)])
            primer_r_seq_rc = reverse_complement(primer_r_seq)
            primer_r_seq_aligned = align_primer(primer_r_seq_rc, amplicon_seq[-len(primer_r_seq_rc):])
            primer_f_mismatch= sum(1 for char in primer_f_seq_aligned if char.islower())
            primer_r_mismatch= sum(1 for char in primer_r_seq_aligned if char.islower())
            original_header = hashmap.get(chrom, '')
            writer.writerow(row + [primer_f_seq, primer_r_seq, amplicon_seq, primer_f_seq_aligned, primer_r_seq_aligned, primer_f_mismatch,primer_r_mismatch, original_header])

def align_primer(primer, amplicon):  # New function**
    aligned_primer = ''
    for p, a in zip(primer, amplicon):
        if match_nucleotides(p, a):
            aligned_primer += p.upper()
        else:
            aligned_primer += p.lower()
            #print("mismatch")
            #print(aligned_primer)
    return aligned_primer

def reverse_complement(seq):  # New function**
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                   'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
                   'H': 'D', 'V': 'B', 'N': 'N'}
    return ''.join(complements[base.upper()] for base in seq[::-1])




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Match primer sequences and amplicon sequences to a table.')
    parser.add_argument('-t', '--table', required=True, help='Path to the input table file.')
    parser.add_argument('-p', '--primers', required=True, help='Path to the primer FASTA file.')
    parser.add_argument('-g', '--genome', required=True, help='Path to the genome FASTA file.')
    parser.add_argument('-m', '--hashmap', required=True, help='Path to the hashmap file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output table file.')
    args = parser.parse_args()

    match_primers_and_amplicons(args.table, args.primers, args.genome,args.hashmap, args.output)
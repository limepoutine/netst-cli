#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
import os
import string
import argparse
import multiprocessing
from multiprocessing import Pool

__description__ = '''Find similar sequences based on haplotype detection output.'''

ARGS = [
    ('-input', dict(metavar='<str>', type=str, help='''input csv file containing sequences with mixed haplotypes.''', required=True)),
    ('-ref_dir', dict(metavar='<str>', type=str, help='''directory containing reference sequences of each haplotype''', required=False, default="")),
    ('-ref_file', dict(metavar='<str>', type=str, help='''file containing reference sequences database''', required=False, default="")),
    ('-k', dict(metavar='<int>', type=int, help='''k-mer size''', required=False, default=21)),
    ('-max', dict(metavar='<int>', type=int, help='''maximum number of similar sequences''', required=False, default=100)),
    ('-out_dir', dict(metavar='<str>', type=str, help='''output directory''', required=True)),
    ('-t', dict(metavar='<int>', type=int, help='''number of concurrent processes''', required=False, default=10))
]

WINDOWS_INVALID_CHARS = '\\/:*?"<>|\''

def generate_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def get_similar_seq(sequence, k, kmer_set):
    count = 0
    for i in range(len(sequence) - k + 1):
        if sequence[i:i+k] in kmer_set:
            count += 1
    return count

def find_matches_from_dir(type_name, k_value, kmer_set, ref_dir):
    file_name = os.path.join(ref_dir, f"{type_name}.fasta")
    sequences_list = []
    with open(file_name, "r") as current_file:
        for line in current_file:
            line = line.strip()
            if line.startswith(">"):
                current_name = line[1:]
            else:
                line = line.upper()
                count = get_similar_seq(line, k_value, kmer_set)
                if count > 0:
                    sequences_list.append([current_name, line, count])
    return sequences_list

def generate_ref_dict(ref_file):
    dic = {}
    file_size = os.path.getsize(ref_file) // 1024
    with open(ref_file, "rb") as current_file:
        while True:
            buf = current_file.peek()
            if not buf:
                break

            empty_chars = 0
            for c in buf:
                if chr(c) not in string.whitespace:
                    break
                empty_chars += 1

            current_file.seek(empty_chars, 1)
            current_pos = current_file.tell()
            line = current_file.readline()
            if not line:
                break
            line = line.decode('utf-8').strip()

            print(f'{current_pos//1024}KB/{file_size}KB', end='\r')

            if line.startswith(">"):
                type_name = line[1:].split("|")[1]
                dic.setdefault(type_name, []).append(current_pos)
    return dic

def find_matches_from_ref_file(type_name, k_value, kmer_set, ref_dict, ref_file):
    ref = ref_dict[type_name]
    sequences_list = []
    with open(ref_file, "r") as current_file:
        for pos in ref:
            current_file.seek(pos)
            current_name = current_file.readline().strip()[1:]
            ref_seq = current_file.readline().strip().upper()
            count = get_similar_seq(ref_seq, k_value, kmer_set)
            if count > 0:
                sequences_list.append([current_name, ref_seq, count])
    return sequences_list

def process_row(row, k_value, out_dir, ref_dir, ref_dict, ref_file, max_sequence):
    seq_id = row[0]
    seq_name = row[1]
    type_list = row[3].split('|')
    support_list = row[4].split('|')
    seq_seq = row[5]
    kmer_set = set(generate_kmers(seq_seq, k_value))
    print(str(seq_id), seq_name, '\t'*8, end='\r')
    with open(os.path.join(out_dir, seq_id + ".fasta"), 'w') as output_file:
        output_file.write(f">{seq_name}|unknown|1\n{seq_seq}\n")
        for i in range(len(type_list)):
            if os.path.isdir(ref_dir):
                sequences_list = find_matches_from_dir(type_list[i], k_value, kmer_set, ref_dir)
            else:
                sequences_list = find_matches_from_ref_file(type_list[i], k_value, kmer_set, ref_dict, ref_file)

            sequences_list = sorted(sequences_list, key=lambda x: x[2], reverse=True)
            for i in range(min(max_sequence,len(sequences_list))):
                output_file.write(f">{sequences_list[i][0]}|{sequences_list[i][2]}\n{sequences_list[i][1]}\n")

def run(k_value, input_file, ref_dir, ref_file, max_sequence, out_dir, num_processes):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    rows = []
    with open(input_file, 'r') as file:
        # Create a CSV reader
        reader = csv.reader(file, delimiter=',')
        next(reader)
        # Iterate over each row in the CSV file
        for row in reader:
            rows.append(row)
            # Print the row

    if os.path.isdir(ref_dir):
        ref_dict = None
    else:
        print("Scanning reference sequences...")
        ref_dict = generate_ref_dict(ref_file)
        print(' '*20, end='\r')

    pool = Pool(processes=num_processes)
    pool.starmap(process_row, [(row, k_value, out_dir, ref_dir, ref_dict, ref_file, max_sequence) for row in rows])
    pool.close()
    pool.join()

def main(pars, args):
    k_value = args.k
    input_file = args.input
    ref_dir = args.ref_dir
    ref_file = args.ref_file
    out_dir = args.out_dir
    max_sequence = args.max
    num_processes = args.t

    if not ref_dir and not ref_file:
        print("Please specify either ref_dir or ref_file.")
        print()
        pars.print_help()
        sys.exit(1)

    run(k_value, input_file, ref_dir, ref_file, max_sequence, out_dir, num_processes)

if __name__ == "__main__":
    if sys.platform.startswith('win'):
        multiprocessing.freeze_support()

    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

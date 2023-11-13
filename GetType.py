#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import csv
import gzip
import argparse

__description__ = '''Detect haplotypes based on a reference database. Mixed haplotypes are detected and support for each putative source haplotype is calculated.'''

ARGS = [
    ('-input', dict(metavar='<str>', type=str, help='''fasta sequences awaiting analysis''', required=True)),
    ('-ref', dict(metavar='<str>', type=str, help='''reference fasta sequences for haplotype detection''', required=True)),
    ('-db', dict(metavar='<str>', type=str, help='''location of the k-mer database (will be created if not existent)''', required=False)),
    ('-k', dict(metavar='<int>', type=int, help='''k-mer size''', required=False, default=21)),
    ('-cut', dict(metavar='<int>', type=int, help='''minimum percentage of identical nucleotides to assume a sequence is not recombinant''', required=False, default=85)),
    ('-output', dict(metavar='<str>', type=str, help='''name prefix of output files''', required=True))
]

def generate_kmers(sequence, k, trans = False):
    kmers = []
    if trans:
        sequence = sequence.translate(str.maketrans('ACGT', '0110'))
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)
    return kmers

def count_kmers_from_file(fasta_file):
    print("Loading k-mer database ...")
    kmer_dict = {}
    with gzip.open(fasta_file, "rt") as file:  # Open the file in gzip mode for reading
        for line in file:
            line = line.strip()
            kmer, sequence_name = line.split("\t")
            kmer_dict[kmer] = sequence_name
    return kmer_dict

def count_kmers(fasta_file, k, output_file):
    if os.path.isfile(output_file):
        kmer_dict  = count_kmers_from_file(output_file)
    else:
        kmer_dict = {}
        current_type = ""
        count = 0
        with open(fasta_file, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    current_type = line[1:].split('|')[1]
                    count += 1
                    print(count, end='\r')
                else:
                    kmers = generate_kmers(line, k)
                    for kmer in kmers:
                        if kmer in kmer_dict:
                            kmer_dict[kmer].add(current_type)
                        else:
                            kmer_dict[kmer] = {current_type}

        for kmer in list(kmer_dict.keys()):
            if len(kmer_dict[kmer]) >= 2:
                kmer_dict.pop(kmer)
            else:
                kmer_dict[kmer] = list(kmer_dict[kmer])[0]
        print("writing kmer dict ...")

        db_dir = os.path.dirname(output_file)
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

        with gzip.open(output_file, "wt") as file:  # Open the file in gzip mode for writing
            for kmer, sequence_name in kmer_dict.items():
                file.write(f"{kmer}\t{sequence_name}\n")

    return kmer_dict

def count_kmer_varieties(target_dict, ref_dict):
    variety_dict = {}
    for kmer in target_dict.keys():
        if kmer in ref_dict:
            my_type = ref_dict[kmer]
            if my_type in variety_dict:
                variety_dict[my_type] += 1
            else:
                variety_dict[my_type] = 1
    return variety_dict,  len(target_dict),

def sort_varieties(variety_dict):
    sorted_varieties = sorted(variety_dict.items(), key=lambda x: x[1], reverse=True)
    return sorted_varieties

def count_kmers_target(fasta_file, ref_dict, result_file, cut_value):
    k = len(next(iter(ref_dict)))
    kmer_dict = {}
    sequence_name = ""
    sequence_seq = ""
    seq_length = 0
    my_id = 0
    with open(fasta_file, "r") as file:
        current_line = ''
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_name:
                    kmers = generate_kmers(current_line, k)
                    sequence_seq = current_line
                    seq_length = len(current_line)
                    for kmer in kmers:
                        if kmer in kmer_dict:
                            kmer_dict[kmer] += 1
                        else:
                            kmer_dict[kmer] = 1

                    kmer_varieties, total_no = count_kmer_varieties(kmer_dict, ref_dict)
                    sorted_varieties = sort_varieties(kmer_varieties)
                    total_count = 0
                    for _, count in sorted_varieties:
                        total_count += count
                    with open(result_file + ".csv", "a", newline='') as rf:
                        writer = csv.writer(rf, delimiter=',')
                        for my_type, count in sorted_varieties:
                            writer.writerow([my_id, sequence_name, seq_length, my_type, 
                                            f"{int(count / total_count * 100)}%",
                                            f"{kmer_varieties[my_type]}",
                                            f"{total_no}",
                                            ])
                    if len(sorted_varieties) > 0:
                        if sorted_varieties[0][1] / total_count * 100 < float(cut_value):
                            with open(result_file + "_mix.csv", "a", newline='') as rf:
                                writer = csv.writer(rf, delimiter=',')
                                my_types = []
                                supports = []
                                for my_type, count in sorted_varieties:
                                    my_types.append(my_type)
                                    supports.append(str(int(count / total_count * 100)))
                                writer.writerow([my_id, sequence_name, seq_length, '|'.join(my_types),
                                '|'.join(supports), sequence_seq])
                    else:
                        with open(result_file + "_null.csv", "a", newline='') as rf:
                            writer = csv.writer(rf, delimiter=',')
                            writer.writerow([my_id, sequence_name, seq_length, sequence_seq])
                    kmer_dict = {}
                sequence_name = line[1:]
                print(sequence_name,end='\r')
                current_line = ''
                my_id += 1
            else:
                current_line += line

    my_id += 1
    kmers = generate_kmers(current_line, k)
    sequence_seq = current_line
    seq_length = len(current_line)
    for kmer in kmers:
        if kmer in kmer_dict:
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer] = 1

    kmer_varieties, total_no = count_kmer_varieties(kmer_dict, ref_dict)
    sorted_varieties = sort_varieties(kmer_varieties)
    total_count = 0
    for _, count in sorted_varieties:
        total_count += count
    with open(result_file + ".csv", "a", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        for my_type, count in sorted_varieties:
            writer.writerow([my_id, sequence_name, seq_length, my_type, 
                            f"{int(count / total_count * 100)}%",
                            f"{kmer_varieties[my_type]}",
                            f"{total_no}",
                            ])

    if len(sorted_varieties) > 0:
        if sorted_varieties[0][1] / total_count * 100 < float(cut_value):
            with open(result_file + "_mix.csv", "a", newline='') as rf:
                writer = csv.writer(rf, delimiter=',')
                my_types = []
                supports = []
                for my_type, count in sorted_varieties:
                    my_types.append(my_type)
                    supports.append(str(int(count / total_count * 100)))
                writer.writerow([my_id, sequence_name, seq_length, '|'.join(my_types),
                '|'.join(supports), sequence_seq])
    else:
        with open(result_file + "_null.csv", "a", newline='') as rf:
            writer = csv.writer(rf, delimiter=',')
            writer.writerow([my_id, sequence_name, seq_length, sequence_seq])

def run(k_value, dict_file, ref_file, target_file, result_file, cut_value):
    kmer_dict_ref = count_kmers(ref_file, k_value, dict_file)

    result_dir = os.path.dirname(result_file)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    with open(result_file + ".csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["ID", "name", "length", "Type", "Support", "hit no", "total no"])

    with open(result_file + "_mix.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["ID", "name", "length", "Type", "Support", "sequence"])

    with open(result_file + "_null.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["ID", "name", "length", "sequence"])

    count_kmers_target(target_file, kmer_dict_ref, result_file, cut_value)

def main(pars, args):
    k_value = args.k
    dict_file = args.db
    ref_file = args.ref
    target_file = args.input
    result_file = args.output
    cut_value = args.cut

    if not dict_file:
        dict_file = ref_file + f'.{k_value}.gz'
        print(f'Generating k-mer database at {dict_file}')

    run(k_value, dict_file, ref_file, target_file, result_file, cut_value)

if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

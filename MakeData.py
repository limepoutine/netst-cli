#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
import os
import multiprocessing
import argparse

__description__ = '''Prepare sequences for haplotype detection and network inference.'''

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Unsupported value encountered.')

ARGS = [
    ('-input', dict(metavar='<str>', type=str, help='''input file''', required=True)),
    ('-file_type', dict(metavar='<str>', type=str, help='''file type (csv or fas)''', required=False, default='fas')),
    ('-index_name', dict(metavar='<int>', type=int, help='''index (0-based) for sequence names''', required=False, default=0)),
    ('-index_type', dict(metavar='<int>', type=int, help='''index (0-based) for sequence types''', required=False, default=1)),
    ('-index_seq', dict(metavar='<int>', type=int, help='''index (0-based) for sequence data (csv only)''', required=False, default=2)),
    ('-out_dir', dict(metavar='<str>', type=str, help='''output directory''', required=True)),
    ('-clean', dict(type=str2bool, nargs='?', const=True, help='''whether to remove duplicate sequences''', default=False)),
    ('-valid', dict(type=str2bool, nargs='?', const=True, help='''whether to filter out invalid sequences''', default=False))
]

def save_sequences_to_files(fasta_file, name_id, type_id, out_dir, do_valid):
    sequences = {}
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        type_name = None
        sequence_id = None
        sequence = ""
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if type_name is not None:
                    # 保存之前的序列
                    if type_name in sequences:
                        sequences[type_name].append([sequence,sequence_id])
                    else:
                        sequences[type_name] = [[sequence,sequence_id]]
                    sequence = ""
                # 提取序列名字
                temp_list = line[1:].split("|")
                sequence_id = temp_list[name_id]
                type_name = temp_list[type_id]
            else:
                # 拼接序列
                sequence += line
        
        # 保存最后一个序列
        if type_name is not None:
            if type_name in sequences:
                sequences[type_name].append([sequence,sequence_id])
            else:
                sequences[type_name] = [[sequence,sequence_id]]
    
    # 将序列保存到对应文件
    for sequence_name, sequence_list in sequences.items():
        file_name = os.path.join(out_dir, f"{sequence_name}.fasta")
        with open(file_name, 'w') as output_file:
            for sequence in sequence_list:
                sequence[0] = sequence[0] .upper()
                if (not do_valid) or is_valid_sequence(sequence[0]):
                    output_file.write(f">{sequence[1]}|{sequence_name}\n")
                    output_file.write(sequence[0]  + "\n")
    return sequences

def is_sequence_contained(sequence, other_sequences):
    for other_seq in other_sequences:
        if other_seq is not sequence and sequence in other_seq:
            return True
    return False

def filter_sequences(input_file, output_file):
    sequences = {}
    with open(input_file, 'r') as f:
        sequence_id = ''
        for line in f:
            if line.startswith('>'):
                sequence_id = line.strip()[1:]
                sequences[sequence_id] = ''
            else:
                sequences[sequence_id] += line.strip()

    filtered_sequences = []
    for seq_id, sequence in sequences.items():
        if not is_sequence_contained(sequence, sequences.values()):
            filtered_sequences.append((seq_id, sequence))

    with open(output_file, 'w') as f:
        for seq_id, sequence in filtered_sequences:
            f.write(f'>{seq_id}\n{sequence}\n')

def process_file(file_path):
    input_file = file_path
    output_file = file_path
    filter_sequences(input_file, output_file)
    print(input_file, "done")

def is_valid_sequence(sequence):
    degenerate_bases = ["N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "-"]
    for base in degenerate_bases:
        if base in sequence:
            return False
    return True

def run(input_file, name_id, type_id, seq_id, out_dir, do_clean, do_valid, file_type):
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    sequences = {}
    if file_type == "csv":
        with open(input_file, 'r', encoding='UTF-8') as file:
            # Create a CSV reader
            reader = csv.reader(file, delimiter=',')
            next(reader)
            # Iterate over each row in the CSV file
            for row in reader:
                # Print the row
                type_name = row[type_id]
                if type_name is not None:
                    # 保存之前的序列
                    if type_name in sequences:
                        sequences[type_name].append([row[seq_id],row[name_id]])
                    else:
                        sequences[type_name] = [[row[seq_id],row[name_id]]]

        # 将序列保存到对应文件

        for type_name, sequence_list in sequences.items():
            print(type_name, end='\r')
            file_name = os.path.join(out_dir, f"{type_name}.fasta")
            with open(file_name, 'w') as output_file:
                for sequence in sequence_list:
                    sequence[0] = sequence[0] .upper()
                    if is_valid_sequence(sequence[0]):
                        output_file.write(f">{sequence[1]}|{type_name}\n")
                        output_file.write(sequence[0]  + "\n")
    else:
        sequences = save_sequences_to_files(input_file, name_id, type_id, out_dir, do_valid)

    if do_clean:
        file_paths = [os.path.join(out_dir, f"{str(i)}.fasta") for i in sequences.keys()]
        pool = multiprocessing.Pool()
        pool.map(process_file, file_paths)
        pool.close()
        pool.join()

    file_paths = [os.path.join(out_dir, f"{str(i)}.fasta") for i in sequences.keys()]
    with open(os.path.join(out_dir, "ref_combine.fasta"), 'w') as combine_file:
        for file_name in file_paths:
            with open(file_name, 'r') as in_file:
                combine_file.write(in_file.read())

def main(pars, args):
    # Extract the file name from the command-line arguments
    input_file = args.input
    name_id = args.index_name
    type_id = args.index_type
    seq_id = args.index_seq
    out_dir = args.out_dir
    do_clean = args.clean
    file_type = args.file_type
    do_valid = args.valid

    run(input_file, name_id, type_id, seq_id, out_dir, do_clean, do_valid, file_type)

if __name__ == "__main__":
    if sys.platform.startswith('win'):
        multiprocessing.freeze_support()

    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

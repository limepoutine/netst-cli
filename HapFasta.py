#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv

__description__ = '''Prepare traits for building haplotype network.'''

ARGS = [
    ('-input', dict(metavar='<str>', type=str, help='''aligned fasta sequences''', required=True)),
    ('-output', dict(metavar='<str>', type=str, help='''name prefix of output files''', required=True))
]

def run(input_file, result_file):
    names = []
    org_seq = []
    AA = "ACGTU"

    new_line = False

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                new_line = True
                names.append(line[1:])
            else:
                if new_line:
                    org_seq.append(line.upper())
                    new_line = False
                else:
                    org_seq[len(org_seq) - 1] += line.upper()

    if not org_seq:
        print("No valid sequence found for building haplotype table.")
        return

    org_seq_clean = org_seq.copy()
    is_var = bytearray(max(map(len, org_seq)))

    for i, seq in enumerate(org_seq):
        print(f"Checking wobbles in seq {names[i]}", end='\r')
        for j, c in enumerate(seq):
            is_var[j] = is_var[j] or (c not in AA)

    for i in range(len(org_seq_clean)):
        temp = list(org_seq_clean[i])
        print(f"Cleaning sequence {names[i]}", end='\r')
        for j in range(len(temp)):
            if is_var[j]:
                temp[j] = '#'
        org_seq_clean[i] = ''.join(temp).replace('#', '')

    hap_seq = []
    hap_name = []

    for i, line in enumerate(org_seq_clean):
        print(f"Calculating haplotype for {names[i]}", end='\r')
        try:
            hap_index = hap_seq.index(line)
            hap_name.append(hap_index)
        except ValueError:
            hap_seq.append(line)
            hap_name.append(len(hap_seq) - 1)

    with open(result_file + "_seq2hap.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["id", "hap", "name", "trait"])
        for i, line in enumerate(org_seq_clean):
            name_trait = names[i].upper().split('|', 1)
            writer.writerow([str(i),
                             f'Hap_{hap_name[i]}',
                             name_trait[0],
                             name_trait[1] if len(name_trait) > 1 else '[NULL]'])

    traits = set()

    for i, name in enumerate(names):
        name_trait = name.upper().split('|', 1)
        s = name_trait[1] if len(name_trait) > 1 else '[NULL]'
        traits.add(s)

    hap_seq = [""] + hap_seq
    trait_list = [""] + list(sorted(traits))

    hap_trait_table = [[0 for i in range(len(trait_list))] for j in range(len(hap_seq))]
    for i in range(1, len(trait_list)):
        hap_trait_table[0][i] = trait_list[i]
    for i in range(1, len(hap_seq)):
        hap_trait_table[i][0] = f'Hap_{hap_name[i - 1]}'
    hap_trait_table[0][0] = ""

    for i, name in enumerate(names):
        name_trait = name.upper().split('|', 1)
        s = name_trait[1] if len(name_trait) > 1 else '[NULL]'
        hap_trait_table[hap_name[i] + 1][trait_list.index(s)] += 1

    with open(result_file + "_hap_trait.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        for row in hap_trait_table:
            writer.writerow(row)

    with open(result_file + "_seq_trait.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(trait_list)
        for name in names:
            name_trait = name.upper().split('|', 1)
            s = name_trait[1] if len(name_trait) > 1 else '[NULL]'
            writer.writerow([name_trait[0]] + [int(t in s.upper()) for t in trait_list[1:]])

    reduced_seq = hap_seq.copy()
    is_var = bytearray(max(map(len, hap_seq[1:])))

    print("Checking haplotype site ...")

    for j in range(len(is_var)):
        is_var[j] = (
            all(len(seq) > j and (seq[j] in AA) for seq in hap_seq[1:]) and
            len(set(seq[j] for seq in hap_seq[1:])) > 1)

    for i in range(1, len(reduced_seq)):
        temp = list(reduced_seq[i])
        for j in range(len(temp)):
            if not is_var[j]:
                temp[j] = '#'
        reduced_seq[i] = ''.join(temp).replace('#', '')

    reduced_seq1 = org_seq_clean.copy()
    for i in range(len(reduced_seq1)):
        temp = list(reduced_seq1[i])
        for j in range(len(temp)):
            if not is_var[j]:
                temp[j] = '#'
        reduced_seq1[i] = ''.join(temp).replace('#', '')

    with open(result_file + "_hap.fasta", "w") as rf:
        for i in range(1, len(reduced_seq)):
            rf.write(f'>Hap_{i - 1}\n')
            rf.write(f'{reduced_seq[i]}\n')

    with open(result_file + "_hap.phy", "w") as rf:
        rf.write(f'{len(reduced_seq) - 1} {len(reduced_seq[1])}\n')
        for i in range(1, len(reduced_seq)):
            rf.write(f'Hap_{i - 1} {reduced_seq[i]}\n')

    with open(result_file + "_seq.phy", "w") as rf:
        rf.write(f'{len(reduced_seq1)} {len(reduced_seq1[0])}\n')
        for i, seq in enumerate(reduced_seq1):
            name = names[i].upper().split('|', 1)[0]
            rf.write(f'{name} {seq}\n')

    with open(result_file + "_seq.fasta", "w") as rf:
        for i, seq in enumerate(reduced_seq1):
            name = names[i].upper().split('|', 1)[0]
            rf.write(f'>{name}\n')
            rf.write(f'{seq}\n')

    with open(result_file + ".meta", "w") as rf:
        for i, seq in enumerate(reduced_seq1):
            name_trait = names[i].upper().split('|', 1)
            s = name_trait[1] if len(name_trait) > 1 else '[NULL]'
            rf.write(f'{name_trait[0]}	{s}\n')

    print(f'Finished creating haplotypes, {len(hap_seq) - 1} haplotypes and {len(reduced_seq[1])} variation sites in total.')

def main(pars, args):
    input_file = args.input
    result_file = args.output

    run(input_file, result_file)

if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

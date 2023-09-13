#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import platform
import shutil
import subprocess
import sys

import GenNetworkConfig
import HapFasta

__description__ = '''Build and visualize a haplotype network from sequences.'''

ARGS = [
    ('-input', dict(metavar='<str>', type=str, help='''sequences to analyze''', required=True)),
    ('-net_type', dict(choices=['original_tcs', 'modified_tcs', 'mjn', 'msn'], help='''type of haplotype network''', required=True)),
    ('-file_type', dict(metavar='<str>', type=str, help='''file type (csv or fas)''', required=False, default='fas')),
    ('-index_name', dict(metavar='<int>', type=int, help='''index (0-based) for sequence names''', required=False, default=1)),
    ('-index_seq', dict(metavar='<int>', type=int, help='''index (0-based) for sequence data (csv only)''', required=False, default=2)),
    ('-index_count', dict(metavar='<int>', type=int, help='''index (0-based) for population count''', required=False, default=-1)),
    ('-index_type', dict(metavar='<int>', type=int, help='''index (0-based) for qualitative (discrete) trait''', required=False, default=-1)),
    ('-index_quant', dict(metavar='<int>', type=int, help='''index (0-based) for quantitative (continuous) trait''', required=False, default=-1)),
    ('-output', dict(metavar='<str>', type=str, help='''output directory''', required=True)),
    ('-aligned', dict(action='store_true', help='''assume input sequences are aligned'''))
]

def csv2fas_generator(csv_file, index_name, index_seq, index_count, index_type, index_quant):
    get_name  = lambda row: row[index_name]
    get_seq   = lambda row: row[index_seq]
    get_type  = lambda row: f'|{row[index_type]}' if index_type != -1 else ''

    if index_count == -1 and index_quant == -1:
        get_traits = lambda row: ''
    else:
        get_count  = lambda row: row[index_count] if index_count != -1 else '1'
        get_quant  = lambda row: row[index_quant] if index_quant != -1 else '1'
        get_traits = lambda row: f'={get_count(row)}={get_quant(row)}'

    with open(csv_file, 'r', encoding='UTF-8') as file:
        reader = csv.reader(file, delimiter=',')
        next(reader)
        for row in reader:
            yield f'>{get_name(row)}{get_traits(row)}{get_type(row)}\n{get_seq(row)}\n'

def fas2fas_generator(fasta_file, index_name, index_seq, index_count, index_type, index_quant):
    get_name  = lambda row: row[index_name]
    get_type  = lambda row: f'|{row[index_type]}' if index_type != -1 else ''

    if index_count == -1 and index_quant == -1:
        get_traits = lambda row: ''
    else:
        get_count  = lambda row: row[index_count] if index_count != -1 else '1'
        get_quant  = lambda row: row[index_quant] if index_quant != -1 else '1'
        get_traits = lambda row: f'={get_count(row)}={get_quant(row)}'

    with open(fasta_file, 'r', encoding='UTF-8') as file:
        row = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if row:
                    yield f'>{get_name(row)}{get_traits(row)}{get_type(row)}\n{sequence}\n'
                    sequence = ""
                row = line[1:].split("|")
            else:
                sequence += line

        yield f'>{get_name(row)}{get_traits(row)}{get_type(row)}\n{sequence}\n'

def run(input_file, aligned, net_type, work_dir):
    script_dir = sys.path[0]

    if not os.path.isdir(script_dir):
        script_dir = os.path.dirname(script_dir)

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    machine = platform.machine()

    if aligned:
        hap_aln_path = input_file
    else:
        hap_aln_path = os.path.join(work_dir, 'hap_aln.fasta')

        if os.path.isfile(hap_aln_path):
            os.remove(hap_aln_path)

        muscle_path = None

        if machine == "x86_64":
            if sys.platform.startswith('win'):
                muscle_path = os.path.join(script_dir, 'muscle5.1.win64.exe')
            elif sys.platform == 'darwin':
                muscle_path = os.path.join(script_dir, 'muscle5.1.macos_intel64')
            else:
                muscle_path = os.path.join(script_dir, 'muscle5.1.linux_intel64')
        elif machine == "arm64" and sys.platform == 'darwin':
            muscle_path = os.path.join(script_dir, 'muscle5.1.macos_arm64')

        if not muscle_path:
            print("Unknown CPU architecture (32-bit is not supported)")
            return

        print("Running muscle 5.1 ...")

        subprocess.run([muscle_path, "-super5", input_file, "-output", hap_aln_path], stdout=sys.stdout, stderr=sys.stderr)

        if not os.path.isfile(hap_aln_path):
            print("Error running muscle 5.1 multiple sequence alignment.")
            return

    net_path = os.path.join(work_dir, 'network')
    seq_phy_path = os.path.join(work_dir, 'network_seq.phy')
    fasthan_out_path = os.path.join(work_dir, 'fasthan')
    gml_path = os.path.join(work_dir, 'fasthan.gml')
    json_path = os.path.join(work_dir, 'fasthan.json')
    meta_path = os.path.join(work_dir, 'network.meta')
    gnn_out_path = os.path.join(work_dir, 'network-config')

    HapFasta.run(hap_aln_path, net_path)

    print("Running FastHaN ...")
    print()

    if os.path.isfile(gml_path):
        os.remove(gml_path)

    fasthan_path = None

    if machine == "x86_64":
        if sys.platform.startswith('win'):
            fasthan_path = os.path.join(script_dir, 'fastHaN_win_intel.exe')
        elif sys.platform == 'darwin':
            pass
        else:
            fasthan_path = os.path.join(script_dir, 'fastHaN_linux')
    elif machine == "arm64":
        if sys.platform.startswith('win'):
            fasthan_path = os.path.join(script_dir, 'fastHaN_win_arm.exe')
        elif sys.platform == 'darwin':
            fasthan_path = os.path.join(script_dir, 'fastHaN_osx_arm')

    if not fasthan_path:
        print("Unknown CPU architecture (32-bit is not supported)")
        return

    subprocess.run([fasthan_path, net_type, "-i", seq_phy_path, "-o", fasthan_out_path], stdout=sys.stdout, stderr=sys.stderr)

    if os.path.isfile(gml_path):
        GenNetworkConfig.run(gml_path, json_path, meta_path, gnn_out_path)

        js_path = os.path.join(work_dir, 'js')
        css_path = os.path.join(work_dir, 'css')
        viewer_path = os.path.abspath(os.path.join(work_dir, 'viewer.html'))

        if not os.path.exists(js_path):
            shutil.copytree(os.path.join(script_dir, 'tcsBU', 'js'), js_path)

        if not os.path.exists(css_path):
            shutil.copytree(os.path.join(script_dir, 'tcsBU', 'css'), css_path)

        if not os.path.exists(viewer_path):
            shutil.copy(os.path.join(script_dir, 'tcsBU', 'index.html'), viewer_path)

        print()
        print(f'Analysis complete. To view the haplotype network, open {viewer_path} in your browser.')
    else:
        print("Error running FastHaN haplotype network analysis.")

def main(pars, args):
    input_file = args.input
    aligned = args.aligned
    net_type = args.net_type
    work_dir = args.output

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    file_type = args.file_type
    index_name = args.index_name
    index_seq = args.index_seq
    index_count = args.index_count
    index_type = args.index_type
    index_quant = args.index_quant

    if file_type == 'csv':
        hap_seq_path = os.path.join(work_dir, 'hap_seq.fasta')

        gen = csv2fas_generator(input_file, index_name, index_seq, index_count, index_type, index_quant)

        with open(hap_seq_path, 'w') as output_file:
            for entry in gen:
                output_file.write(entry)

        input_file = hap_seq_path
    elif index_count != -1 or index_type != -1 or index_quant != -1:
        hap_seq_path = os.path.join(work_dir, 'hap_seq.fasta')

        gen = fas2fas_generator(input_file, index_name, index_seq, index_count, index_type, index_quant)

        with open(hap_seq_path, 'w') as output_file:
            for entry in gen:
                output_file.write(entry)

        input_file = hap_seq_path

    run(input_file, aligned, net_type, work_dir)

if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import datetime
import importlib
import multiprocessing
import os
import sys

__description__ = '''K-mer based haplotype detection and haplotype network building tool.'''

SUBCOMMAND_DESC = f'''These subcommands invoke a specific stage of processing. Use "pipeline" to run the whole pipeline.

makedata    convert a packed csv or fasta file into individual fasta sequences
getblast    infer haplotypes from sequences using blastn search
gettype     infer mixed haplotypes from sequences using k-mer based search
getmixseq   search for fragments supporting each haplotype using gettype results
viewhaplo   build and visualize a haplotype network from sequences
pipeline    run the haplotype network analysis pipeline
'''

MODULE_NAMES = ["MakeData", "GetBlast", "GetType", "GetMixSeq", "ViewHaplo"] # "GenNetworkConfig", "HapFasta"

def run_pipeline(args):
    import GetBlast
    import GetType
    import GetMixSeq
    import GenNetworkConfig
    import ViewHaplo

    work_dir = args.work_dir
    method_type = args.method
    net_type = args.net_type

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    target_file = args.input
    ref_file = args.ref
    cut_value = args.cut

    k_value = args.k

    e_value = args.e
    word_size = args.word_size

    haplo_seq = args.haplo_seq
    haplo_seq_type = args.haplo_seq_type
    index_name = args.index_name
    index_seq = args.index_seq
    index_count = args.index_count
    index_type = args.index_type
    index_quant = args.index_quant

    skip_getmixseq = args.skip_getmixseq
    include_all_seqs = args.all_seqs
    process_num = args.t

    if args.db:
        db_file = args.db
    elif method_type == "blast":
        db_file = os.path.join(work_dir, 'blast', 'db')
    else:
        db_file = os.path.join(work_dir, f'db_{k_value}.gz')

    result_file = os.path.join(work_dir, 'type_analysis')
    mix_csv_path = os.path.join(work_dir, 'type_analysis_mix.csv')
    null_csv_path = os.path.join(work_dir, 'type_analysis_null.csv')
    all_csv_path = os.path.join(work_dir, 'type_analysis.csv')

    if method_type == "blast":
        GetBlast.run(target_file, ref_file, db_file, None, e_value, word_size, result_file, cut_value, process_num)
    else:
        GetType.run(k_value, db_file, ref_file, target_file, result_file, cut_value)

    print("Finished haplotype analysis.")
    print()

    if not skip_getmixseq:
        current_time = datetime.datetime.now().strftime("%y%m%d%H%M%S")

        GetMixSeq.run(k_value,
                      mix_csv_path,
                      "",
                      ref_file,
                      100,
                      os.path.join(work_dir, f'mix-{current_time}'),
                      process_num)

        print("Finished sequence search.")
        print()

    hap_seq_path = os.path.join(work_dir, 'hap_seq.fasta')
    network_dir = os.path.join(work_dir, 'network')

    single_haps = {}
    mix_haps = set()

    with open(hap_seq_path, 'w') as wfile:
        if haplo_seq:
            if haplo_seq_type == 'csv':
                c2f = ViewHaplo.csv2fas_generator(haplo_seq, index_name, index_seq, index_count, index_type, index_quant)
                for line in c2f:
                    wfile.write(line)
            elif index_count != -1 or index_type != -1 or index_quant != -1:
                f2f = ViewHaplo.fas2fas_generator(haplo_seq, index_name, index_seq, index_count, index_type, index_quant)
                for line in f2f:
                    wfile.write(line)
            else:
                with open(haplo_seq, 'r') as rfile:
                    for line in rfile:
                        wfile.write(line)

        with open(mix_csv_path, 'r') as rfile:
            reader = csv.reader(rfile, delimiter=',')
            next(reader)
            for row in reader:
                mix_haps.add(row[1])
                wfile.write(f">{row[1]}|{row[3]}\n")
                wfile.write(f"{row[5]}\n")

        with open(null_csv_path, 'r') as rfile:
            reader = csv.reader(rfile, delimiter=',')
            next(reader)
            for row in reader:
                mix_haps.add(row[1])
                wfile.write(f">{row[1]}|[null]\n")
                wfile.write(f"{row[3]}\n")

        if include_all_seqs:
            with open(all_csv_path, 'r') as rfile:
                reader = csv.reader(rfile, delimiter=',')
                next(reader)
                for row in reader:
                    if row[3] in single_haps:
                        mix_haps.add(row[1])
                    else:
                        single_haps[row[1]] = row[3]

            for name in mix_haps:
                if name in single_haps:
                    del single_haps[name]

            seq_name = ''

            with open(target_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        seq_name = line[1:]
                        if seq_name in single_haps:
                            wfile.write(f">{seq_name}|{single_haps[seq_name]}\n")
                    elif seq_name in single_haps:
                        wfile.write(f"{line}\n")

    if not single_haps and not mix_haps:
        print("No haplotypes eligible for haplotype network. You might want to run with -all_seqs to include all sequences in the network.")
        return

    ViewHaplo.run(hap_seq_path, False, net_type, network_dir)

if __name__ == "__main__":
    if sys.platform.startswith('win'):
        multiprocessing.freeze_support()

    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    subpars = pars.add_subparsers(dest='module_name', description=SUBCOMMAND_DESC, metavar='<command>', help='one of the above commands')

    sub_parsers = {}
    modules = {}

    # Load each script into modules dict (with keys in lower case)
    for name in MODULE_NAMES:
        mod = importlib.import_module(name)
        mod_pars = subpars.add_parser(name.lower(),
                                      formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description=mod.__description__)

        for param in mod.ARGS:
            mod_pars.add_argument(param[0], **param[1])

        modules[name.lower()] = mod
        sub_parsers[name.lower()] = mod_pars

    pip_pars = subpars.add_parser('pipeline', description='Run the haplotype analysis pipeline. A haplotype network will be generated from analyzed sequences and the reference database.')

    pip_pars.add_argument('-work_dir', metavar='<str>', type=str, help='''directory to save intermediate files and analysis results''', required=True)
    pip_pars.add_argument('-method', choices=['blast', 'kmer'], help='''haplotype detection method (default kmer)''', required=False, default='kmer')
    pip_pars.add_argument('-net_type', choices=['original_tcs', 'modified_tcs', 'mjn', 'msn'], help='''type of haplotype network (default modified_tcs)''', required=False, default='modified_tcs')

    pip_pars.add_argument('-input', metavar='<str>', type=str, help='''fasta sequences awaiting analysis''', required=True)
    pip_pars.add_argument('-ref', metavar='<str>', type=str, help='''reference fasta sequences for haplotype detection''', required=True)
    pip_pars.add_argument('-db', metavar='<str>', type=str, help='''location of pre-built database for the selected detection method''', required=False)
    pip_pars.add_argument('-cut', metavar='<int>', type=int, help='''minimum percentage of identical nucleotides to assume a sequence is not recombinant''', required=False, default=85)

    pip_pars.add_argument('-k', metavar='<int>', type=int, help='''k-mer size''', required=False, default=21)

    pip_pars.add_argument('-e', metavar='<num>', type=float, help='''maximum E value in BLAST queries''', required=False, default=1e-5)
    pip_pars.add_argument('-word_size', metavar='<int>', type=int, help='''word size in BLAST queries''', required=False, default=11)

    pip_pars.add_argument('-haplo_seq', metavar='<str>', type=str, help='''other sequences to display in haplotype network''', required=False)
    pip_pars.add_argument('-haplo_seq_type', metavar='<str>', type=str, help='''type of haplo_seq (csv or fas)''', required=False, default='fas')
    pip_pars.add_argument('-index_name', metavar='<int>', type=int, help='''index (0-based) for names in haplo_seq''', required=False, default=1)
    pip_pars.add_argument('-index_seq', metavar='<int>', type=int, help='''index (0-based) for sequence data in haplo_seq (csv only)''', required=False, default=2)
    pip_pars.add_argument('-index_count', metavar='<int>', type=int, help='''index (0-based) for population count in haplo_seq''', required=False, default=-1)
    pip_pars.add_argument('-index_type', metavar='<int>', type=int, help='''index (0-based) for haplotype in haplo_seq''', required=False, default=-1)
    pip_pars.add_argument('-index_quant', metavar='<int>', type=int, help='''index (0-based) for quantitative (continuous) trait''', required=False, default=-1)

    pip_pars.add_argument('-all_seqs', help='''include putative non-recombinant sequences in haplotype network''', action='store_true')
    pip_pars.add_argument('-skip_getmixseq', help='''skip time-consuming search of haplotype-supporting fragments''', action='store_true')
    pip_pars.add_argument('-t', metavar='<int>', type=int, help='''number of concurrent processes''', required=False, default=10)

    args = pars.parse_args()

    if not args.module_name:
        # No subcommand, print help and exit
        print("Please specify a subcommand to execute.")
        print()
        pars.print_help()
        sys.exit(1)

    module_name = args.module_name.lower()

    if module_name != "pipeline":
        # Invoke the respective module with args
        modules[module_name].main(sub_parsers[module_name], args)
    else:
        run_pipeline(args)

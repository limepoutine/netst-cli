#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import collections
import csv
import os
import string
import subprocess
import sys

__description__ = '''Detect haplotypes using BLAST queries.'''

ARGS = [
    ('-input', dict(metavar='<str>', type=str, help='''fasta sequences awaiting analysis''', required=True)),
    ('-ref', dict(metavar='<str>', type=str, help='''reference fasta sequences for haplotype detection''', required=True)),
    ('-db', dict(metavar='<str>', type=str, help='''location of the BLAST database (will be created if not existent)''', required=False)),
    ('-blast_result', dict(metavar='<str>', type=str, help='''location of pre-computed BLAST results (requires -outfmt 6)''', required=False)),
    ('-e', dict(metavar='<num>', type=float, help='''maximum E value in BLAST queries''', required=False, default=1e-5)),
    ('-word_size', dict(metavar='<int>', type=int, help='''word size in BLAST queries''', required=False, default=11)),
    ('-cut', dict(metavar='<int>', type=int, help='''minimum percentage of identical nucleotides to assume a sequence is not recombinant''', required=False, default=85)),
    ('-output', dict(metavar='<str>', type=str, help='''name prefix of output files''', required=True)),
    ('-t', dict(metavar='<int>', type=int, help='''number of concurrent threads''', required=False, default=10))
]

# blastn outfmt 6
# 0     	1     	2     	3     	4       	5      	6     	7   	8     	9   	10    	11
# qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore
BlastResultLine = collections.namedtuple("ResultLine", ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

def get_script_dir():
    script_dir = sys.path[0]

    if not os.path.isdir(script_dir):
        script_dir = os.path.dirname(script_dir)

    return script_dir

def generate_seq_dict(ref_file):
    id_dict = {}
    ptr_dict = {}

    current_id = 0
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
                current_id += 1
                type_name = line[1:]
                id_dict[type_name] = current_id
                ptr_dict[type_name] = current_pos

    return id_dict, ptr_dict

def run_blast(target_file, db_file, blast_result, e_value, word_size, result_file, cut_value, num_threads):
    script_dir = get_script_dir()

    id_dict, seq_dict = generate_seq_dict(target_file)

    if blast_result:
        out_path = blast_result
    else:
        print('Running blastn query ...')

        out_path = result_file + "_blast.tbl"

        if sys.platform.startswith('win'):
            blastn_path = os.path.join(script_dir, 'blastn.exe')
        elif sys.platform == 'darwin':
            raise RuntimeError("TODO") # TODO
        else:
            blastn_path = os.path.join(script_dir, 'blastn')

        proc = subprocess.run([blastn_path,
                               "-query", target_file,
                               "-db", db_file,
                               "-out", out_path,
                               "-task", "dc-megablast",
                               "-word_size", str(word_size),
                               "-evalue", str(e_value),
                               "-outfmt", "6",
                               "-mt_mode", "1",
                               "-num_threads", str(num_threads)],
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)

    results = {}

    with open(out_path, 'r') as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            res = BlastResultLine(*line.split("\t"))
            res = res._replace(
                pident=float(res.pident),
                length=int(res.length),
                mismatch=int(res.mismatch),
                gapopen=int(res.gapopen),
                qstart=int(res.qstart),
                qend=int(res.qend),
                sstart=int(res.sstart),
                send=int(res.send),
                evalue=float(res.evalue),
                bitscore=float(res.bitscore))

            results.setdefault(res.qseqid, []).append(res)

    with open(target_file, "r") as f:
        for seq_name, matches in results.items():
            seq_id = id_dict[seq_name]

            f.seek(seq_dict[seq_name])
            f.readline()

            sequence_seq = ''

            for line in f:
                line = line.strip()

                if line.startswith(">"):
                    break
                else:
                    sequence_seq += line.upper()

            seq_length = len(sequence_seq)

            highest_hap = {}

            for ent in matches:
                name, hap = ent.sseqid.split('|', )

                if hap not in highest_hap:
                    highest_hap[hap] = {'name': name, 'overlap': (ent.qend - ent.qstart) * ent.pident / 100}
                    last_name = name
                elif name == highest_hap[hap]['name']:
                    highest_hap[hap]['overlap'] += (ent.qend - ent.qstart) * ent.pident / 100

            highest_hap = list(sorted(highest_hap.items(), key=lambda x: x[1]['overlap'], reverse=True))
            overlap_sum = sum(map(lambda x: x[1]['overlap'], highest_hap))
            types, supports = zip(*map(lambda x: (x[0], x[1]['overlap'] / overlap_sum), highest_hap))

            with open(result_file + ".csv", "a", newline='') as rf:
                writer = csv.writer(rf, delimiter=',')
                for hap, ent in highest_hap:
                    writer.writerow([seq_id, seq_name, seq_length, hap, int(ent['overlap']), seq_length])

            if supports[0] * 100 < float(cut_value):
                with open(result_file + "_mix.csv", "a", newline='') as rf:
                    writer = csv.writer(rf, delimiter=',')
                    writer.writerow([seq_id, seq_name, seq_length, '|'.join(types),
                                        '|'.join(str(int(s * 100)) for s in supports),
                                        sequence_seq])

        # Those not in blast results
        with open(result_file + "_null.csv", "a", newline='') as rf:
            writer = csv.writer(rf, delimiter=',')
            for seq_name in id_dict:
                seq_id = id_dict[seq_name]

                f.seek(seq_dict[seq_name])
                f.readline()
                sequence_seq = f.readline().strip().upper()
                seq_length = len(sequence_seq)

                if seq_name not in results:
                    writer.writerow([seq_id, seq_name, seq_length, sequence_seq])

def run_makeblastdb(ref_file, db_file):
    script_dir = get_script_dir()

    print('Running makeblastdb ...')

    if sys.platform.startswith('win'):
        makeblastdb_path = os.path.join(script_dir, 'makeblastdb.exe')
    elif sys.platform == 'darwin':
        raise RuntimeError("TODO") # TODO
    else:
        makeblastdb_path = os.path.join(script_dir, 'makeblastdb')

    proc = subprocess.run([makeblastdb_path,
                           "-dbtype", "nucl",
                           "-in", ref_file,
                           "-input_type", "fasta",
                           "-parse_seqids",
                           "-hash_index",
                           "-out", db_file],
                           stdout=sys.stdout, stderr=sys.stderr)

    return proc.returncode

def run(target_file, ref_file, db_file, blast_result, e_value, word_size, result_file, cut_value, num_threads):
    with open(result_file + ".csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["ID", "name", "length", "Type", "Support", "hit len", "total len"])

    with open(result_file + "_mix.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["ID", "name", "length", "Type", "Support", "sequence"])

    with open(result_file + "_null.csv", "w", newline='') as rf:
        writer = csv.writer(rf, delimiter=',')
        writer.writerow(["ID", "name", "length", "sequence"])

    if not blast_result and not os.path.isfile(db_file + ".nhr"):
        ret = run_makeblastdb(ref_file, db_file)

        if ret != 0:
            print("Unable to build BLAST database.")
            return

    run_blast(target_file, db_file, blast_result, e_value, word_size, result_file, cut_value, num_threads)

def main(pars, args):
    target_file = args.input
    ref_file = args.ref
    db_file = args.db
    blast_result = args.blast_result
    result_file = args.output
    e_value = args.e
    word_size = args.word_size
    cut_value = args.cut
    num_threads = args.t

    if not db_file and not blast_result:
        db_file = ref_file + '.blast_db'
        print(f'Making BLAST database at {db_file}')

    run(target_file, ref_file, db_file, blast_result, e_value, word_size, result_file, cut_value, num_threads)

if __name__ == "__main__":
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

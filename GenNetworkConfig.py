    #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#*************************************************************************
#    > File Name: GenNetworkConfig.py
#    > Author: xlzh
#    > Mail: xiaolongzhang2015@163.com 
#    > Created Time: 2021年12月02日 星期四 22时20分30秒
#*************************************************************************
# Modified By YY 2023/7/9

import argparse
import sys
import json
import random

__description__ = '''Generates Javascript for haplotype network viewer.'''

ARGS = [
    ('gml_input', dict(type=str, help='''input gml file''')),
    ('json_input', dict(type=str, help='''input json file''')),
    ('meta_input', dict(type=str, help='''input meta file''')),
    ('out_prefix', dict(type=str, help='''output file name prefix'''))
]

def json_parse(json_file):
    ''' func: parse the json file to obtain the haplotye
    '''
    hap_list = []
    json_wraper = json.load(open(json_file, 'r'))

    # node -> {'id': 0, 'frequency': 1.0, 'title1': 'Sample1', 'title2': 'Sample1;Sample2'}
    for node in json_wraper['nodes']: 
        hap_list.extend([_.split("=")[0] for _ in node['title2'].split(';')])

    return hap_list


def read_meta_file(meta_file):
    ''' func: read the meta file and obtain the country for each individual
        meta_dict = {'EPI_ISL_6814923': 'Australia', ...}
    '''
    meta_dict = {}
    meta_fp = open(meta_file, 'r')

    for line in meta_fp:
        l = line.rstrip().split('\t')
        meta_dict[l[0].split('=')[0]] = l[1]  # l[0] -> sampleid; l[1] -> country

    return meta_dict


def gen_hap_config(hap_list, meta_dict):
    ''' func: generate the haplotype config
        hap_conf_list = [(sample1, China), (sample2, Korean), ...]
    '''
    hap_conf_list = []

    for hap in hap_list:  # EPI_ISL_6832737 or 'IN12'
        if hap.startswith('IN'):  # skip the intermediate node
            continue

        elif hap not in meta_dict:
            sys.stderr.write("[warning] the group of %s is not given in the meta file\n" % hap)

        else:
            country = meta_dict[hap].replace(' ', '')
            hap_conf_list.append((hap, country))

    return hap_conf_list


def _random_color(seed):
    ''' func: generate color randomly based on the specified seed
    '''
    random.seed(seed)

    colorArr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    color_list = []

    for i in range(6):
        color_list.append(colorArr[random.randint(0, 14)])

    return '#' + ''.join(color_list)

def text_to_integer(text):
    # 将文本中的每个字符映射成对应的ASCII码，并拼接在一起
    integer_value = int(''.join(str(ord(char)) for char in text))
    return integer_value

def gen_group_config(hap_conf_list):
    ''' func: generate the group config
        group_conf_list = [(China, '#FF011B', 'none'), (Korean, '2693FF', 'none'), ...]
    '''
    group_conf_list = []
    count_dict = {}

    for item in hap_conf_list:
        if item[1] != 'Default' and item[1] not in count_dict:
            rand_color = _random_color(text_to_integer(item[1]))
            print(text_to_integer(item[1]), rand_color)
            count_dict[item[1]] = rand_color
            group_conf_list.append((item[1], rand_color, 'none'))

    return group_conf_list


def write_conf(hap_conf_list, group_conf_list, out_prefix):
    ''' func: write the config file
    '''
    # write the haplotype config file
    hap_fp = open(out_prefix + '_hapconf.csv', 'w')
    for hap in hap_conf_list:
        hap_fp.write("%s;%s\n" % (hap[0], hap[1]))

    # write the group config file
    group_fp = open(out_prefix + '_groupconf.csv', 'w')
    for group in group_conf_list:
        group_fp.write("%s;%s;%s\n" % (group[0], group[1], group[2]))

    hap_fp.close()
    group_fp.close()


def file2line(file_name):
    content = ""
    with open(file_name, 'r') as file_fp: 
        content = file_fp.read()
        content = content.replace('\"',"\\\"").replace('\n',"\\n")
    return content

def run(gml_file, json_file, meta_file, out_prefix):
    hap_list = json_parse(json_file)
    meta_dict = read_meta_file(meta_file)
    hap_conf_list = gen_hap_config(hap_list, meta_dict)
    group_conf_list = gen_group_config(hap_conf_list)
    write_conf(hap_conf_list, group_conf_list, out_prefix)

    data_fp = open(out_prefix + '.js', 'w')
    data_fp.write("var gmlfile = {target: {files: [new File([\"")
    data_fp.write(file2line(gml_file))
    data_fp.write("\"], \".gml\")]}};\n")

    data_fp.write("var hapconffile = {target: {files: [new File([\"")
    data_fp.write(file2line(out_prefix + '_hapconf.csv'))
    data_fp.write("\"], \".gml\")]}};\n")

    data_fp.write("var groupconffile = {target: {files: [new File([\"")
    data_fp.write(file2line(out_prefix + '_groupconf.csv'))
    data_fp.write("\"], \".gml\")]}};\n")
    data_fp.close()

def main(pars, args):
    gml_file = args.gml_input
    json_file = args.json_input
    meta_file = args.meta_input
    out_prefix = args.out_prefix

    run(gml_file, json_file, meta_file, out_prefix)

if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__description__)
    for param in ARGS:
        pars.add_argument(param[0], **param[1])
    args = pars.parse_args()

    main(pars, args)

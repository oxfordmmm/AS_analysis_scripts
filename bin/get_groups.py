#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import gzip
import toml
from argparse import ArgumentParser
from Bio import SeqIO

def channelsToml(tomlFILE):
    toml_d = toml.load(tomlFILE)
    channels={}
    for condition in toml_d['conditions']:
        for channel in toml_d['conditions'][condition]['channels']:
            channels[channel]=toml_d['conditions'][condition]['name']
        #channels[toml_d['conditions'][condition]['name']]=toml_d['conditions'][condition]['channels']
    
    return channels

def _get_read_info(fq):
    for record in SeqIO.parse(gzip.open(fq, 'rt'), "fastq"):
        desc=record.description.split(' ')
        desc_dict={}
        for d in desc:
            if '=' in d:
                key,value=d.split('=')
                desc_dict[key]=value
        yield record.id, int(desc_dict['ch']), desc_dict['start_time']


def run(opts):
    chs=channelsToml(opts.toml)
    unblocked=open(opts.unblocks).read().split('\n')
    
    AS_unblocked, AS_sequenced, control_sequenced = [],[],[]
    for i in _get_read_info(opts.fastq):
        if i[0] in unblocked:
            if chs[i[1]] != "control condtion 1":
                AS_unblocked.append(i)
            else:
                control_sequenced.append(i)
        else:
            if chs[i[1]] != "control condtion 1":
                AS_sequenced.append(i)
            else:
                control_sequenced.append(i)

    with open('AS_unblocked.txt', 'wt') as f:
        for i in AS_unblocked:
            f.write(f'{i[0]}\n')
    
    with open('AS_sequenced.txt', 'wt') as f:
        for i in AS_sequenced:
            f.write(f'{i[0]}\n')
    
    with open('control_sequenced.txt', 'wt') as f:
        for i in control_sequenced:
            f.write(f'{i[0]}\n')
    

def args(parser):
    parser.add_argument('-fq', '--fastq', required=False,
                             help='fastq file')
    parser.add_argument('-u', '--unblocks', required=False,
                             help='unblocks file')
    parser.add_argument('-c', '--toml', required=True,
                             help='chanels toml file')

    return parser

if __name__=="__main__":
    # args
    parser = ArgumentParser(description='Conditions, classifications and times from centrifuge and seqsums')
    parser = args(parser)
    opts, unknown_args = parser.parse_known_args()
    # run script

    run(opts)

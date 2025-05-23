#!/usr/bin/env python3
from ete3 import NCBITaxa
import sys
import argparse

def run(inputs, output):
    ncbi = NCBITaxa()

    # get descendent taxa for Clostridium difficile
    taxids=ncbi.get_descendant_taxa('186801', intermediate_nodes=True)    

    taxids.append(1496)

    reads=[]
    for i in inputs:
        for line in open(i,'rt'):
            line=line.split('\t')
            if int(line[2]) in taxids:
                reads.append(line[1])

    with open(output,'wt') as f:
        f.write('\n'.join(reads))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get read ids from taxids')
    parser.add_argument('-i','--inputs', nargs='+', help='Input files')
    parser.add_argument('-o', '--output', help='Output file')
    args = parser.parse_args()
    run(args.inputs, args.output)

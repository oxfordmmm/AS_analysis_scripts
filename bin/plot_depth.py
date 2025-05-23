#!/usr/bin/env python3
import pandas as pd
import sys
from argparse import ArgumentParser
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="whitegrid")


def getData(f:str):
    names=['chrom','position','depth']
    df=pd.read_csv(f,sep='\t',names=names)
    df=df.sample(1000)
    return df

def plot_depth(df, output, barcode=None, taxid=None, batch=None):
    '''Draw a line plot of depth with position on the X axis and depth on the Y axis'''
    sns.lineplot(x='position',y='depth',data=df)
    plt.title(f'{barcode} {taxid} {batch}')
    plt.savefig(output)
    
    
if __name__ == '__main__':
    parser = ArgumentParser(description='')
    parser.add_argument('-i', '--input', required=True, 
                             help='Input tsv file from samtools depth')
    parser.add_argument('-o', '--output', required=True, 
                             help='Output csv file summary')
    parser.add_argument('-b', '--barcode', required=True, 
                             help='Barcode name, number')
    parser.add_argument('-a', '--batch', required=True, 
                             help='batch name, number')
    parser.add_argument('-t', '--taxid', required=True, 
                             help='taxid number')
    opts, unknown_args = parser.parse_known_args()

    df=getData(opts.input)
    df=plot_depth(df, opts.output, barcode=opts.barcode, taxid=opts.taxid, batch=opts.batch)

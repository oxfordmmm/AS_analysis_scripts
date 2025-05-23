#!/usr/bin/env python3
import pandas as pd
import sys
from argparse import ArgumentParser
from Bio import SeqIO
import gzip

def cluster_regions(r):
    n=0
    r2=r
    if r==r2+1:
        r2=r
        yield n
    else:
        r2=r
        n+=1
        yield n


def getData(f:str):
    names=['chrom','position','depth']
    df=pd.read_csv(f,sep='\t',names=names)
    df=df[df['depth']>0]

    # count consecutive positions
    clusters=[]
    n=0
    i2=0
    for i in df['position'].to_list():
        if i==i2+1:
            clusters.append(n)
            i2=i
        else:
            n+=1
            i2=i
            clusters.append(n)
    df['cluster']=list(clusters)
    df['cluster length']=df.groupby('cluster')['position'].transform(len)
    df['cluster start']=df.groupby('cluster')['position'].transform('first')
    df['cluster end']=df.groupby('cluster')['position'].transform('last')
    return df


def get_ref_regions(df, ref):
    df=df[['chrom','cluster','cluster start','cluster end']]
    df.drop_duplicates(inplace=True)
    for seq in SeqIO.parse(gzip.open(ref, 'rt'), "fasta"):
        for index, row in df.iterrows():
            if row['chrom']==seq.id:
                print(f'>{seq.id}:{row["cluster start"]}-{row["cluster end"]}')
                print(seq.seq[row['cluster start']:row['cluster end']])


    
if __name__ == '__main__':
    parser = ArgumentParser(description='')
    parser.add_argument('-i', '--input', required=True, 
                             help='Input tsv file from samtools depth')
    parser.add_argument('-r', '--ref', required=True, 
                             help='ref genome')
    parser.add_argument('-n', '--num_consecutive', required=True, 
                             help='number of consecutive positions to consider before extracting region')
    opts, unknown_args = parser.parse_known_args()

    df=getData(opts.input)
    df.to_csv('test.csv',index=False)
    df=df[df['cluster length']>=int(opts.num_consecutive)]
    get_ref_regions(df,opts.ref)
    #df=plot_depth(df, opts.output, barcode=opts.barcode, taxid=opts.taxid, batch=opts.batch)

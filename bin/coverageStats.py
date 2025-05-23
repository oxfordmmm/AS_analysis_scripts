#!/usr/bin/env python3
import pandas as pd
import sys
from argparse import ArgumentParser


def getData(f:str):
    names=['chrom','position','depth']
    df=pd.read_csv(f,sep='\t',names=names)
    return df

def coverageStats(df, barcode=None, taxid=None, batch=None):
    cov1=df[df.depth > 0].groupby(['chrom']).count()
    cov10=df[df.depth > 9].groupby(['chrom']).count()

    bases=df.groupby(['chrom'])['depth'].sum()
    chromLens=df.groupby(['chrom'])['position'].count()
    avDepth=bases/cov1['depth']

    cov1.reset_index(inplace=True)
    cov10.reset_index(inplace=True)

    df2=cov1.merge(cov10,on='chrom',suffixes=[' cov1',' cov10'],how='outer')
    df2['length']=df2.chrom.map(chromLens)
    df2['covBreadth1x']=df2['depth cov1']/df2['length']
    df2['covBreadth10x']=df2['depth cov10']/df2['length']
    df2['avDepth']=df2.chrom.map(avDepth)
    df2['Barcode']=barcode
    df2['taxid']=taxid
    df2['batch']=batch
    df2['bases']=df2.chrom.map(bases)
    df2=df2[['Barcode','batch','taxid', 'chrom','length','bases','avDepth','position cov1','position cov10','covBreadth1x','covBreadth10x']]
    
    #print(df2)
    return df2
    
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
    df=coverageStats(df, barcode=opts.barcode, taxid=opts.taxid, batch=opts.batch)
    df.to_csv(opts.output,index=False)

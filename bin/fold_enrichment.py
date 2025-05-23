#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")

def getFoldEnrichmentTime(f, species, AS_type):
    df=pd.read_csv(f,sep='\t')
    cols=['Condition','Read type','5mins','start_time', 'Group reads per 5mins', 'Condtion reads per 5mins', 'Group cumulative bases','sequence_length_template']
    df2=df[cols]
    df2['Culmunitive max bases 5 mins']=df2.groupby(['Condition','Read type','5mins'])[['Group cumulative bases']].transform(max)
    df2['Culmunitive median bases 5 mins']=df2.groupby(['Condition','Read type','5mins'])[['Group cumulative bases']].transform(np.median)

    df3=df2.drop_duplicates(subset=['Condition', 'Read type', '5mins'])
    dfAS=df3[(df3['Condition']==AS_type) & (df3['Read type']==species)]
    dfControl=df3[(df3['Condition']=='control condtion 1') & (df3['Read type']==species)]
    df=dfAS.merge(dfControl, on='5mins', how='outer', suffixes=('_AS','_control'))
    df.fillna(method='ffill', inplace=True)

    df['Fold enrichment']=df['Culmunitive max bases 5 mins_AS']/df['Culmunitive max bases 5 mins_control']
    df['Fold enrichment median 5 mins']=df['Culmunitive median bases 5 mins_AS']/df['Culmunitive median bases 5 mins_control']
    df['Fold enrichment raw']=df['Group cumulative bases_AS']/df['Group cumulative bases_control']


    return df




runs={"AS_S1_Ecoli_10_4_deplete_0": {'species': 'E.coli', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'Run 1', 'library':1},
"AS_S1_Ecoli_10_4_enrich_0": {'species': 'E.coli', 'AS_type': 'Enrich E.coli', 'Actual AS type': 'Enrich', 'Alt run name':'Run 2', 'library':1},
"AS_S2_1_Ecol_10_4_deplete_0": {'species': 'E.coli', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'Run 3', 'library':2},
"AS_S2_1_Ecol_10_4_enrich_0": {'species': 'E.coli', 'AS_type': 'Enrich E.coli', 'Actual AS type': 'Enrich', 'Alt run name':'Run 4', 'library':2},
"AS_S2_1_Saur_10_4_enrich_0": {'species': 'S.aureus', 'AS_type': 'Enrich E.coli', 'Actual AS type': 'Enrich', 'Alt run name':'Run 5', 'library':3},
"AS_S2_2_Ecol_10_4_deplete_0": {'species': 'E.coli', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'Run 6', 'library':4},
"AS_S2_2_Ecol_10_4_enrich_0": {'species': 'E.coli', 'AS_type': 'Enrich E.coli', 'Actual AS type': 'Enrich', 'Alt run name':'Run 7', 'library':4},
"AS_S2_2_Saur_10_4_enrich_0": {'species': 'S.aureus', 'AS_type': 'Enrich E.coli', 'Actual AS type': 'Enrich', 'Alt run name':'Run 8', 'library':5},
"AS_Staph_50_50_deplete_0": {'species': 'S.aureus', 'AS_type': 'Deplete Human', 'Actual AS type': 'Deplete', 'Alt run name':'Run 9', 'library':6},
"AS_T1_deplete_0": {'species': 'E.faecalis', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'Run 10', 'library':7},
"saur50_human50_20k_ligation": {'species': 'S.aureus', 'AS_type': 'enrich target', 'Actual AS type': 'Deplete', 'Alt run name':'Run 11', 'library':8}}

runs={"sample_1a": {'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_1a', 'library':9},
"sample_2a": {'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_2a', 'library':9},
"sample_3a": {'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_3a', 'library':9},
"sample_4a": {'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_4a', 'library':9},
"sample_6a":{'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_6a', 'library':9},
"sample_7a":{'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_7a', 'library':9},
"sample_PBS_neg1": {'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_PBS_neg1', 'library':9},
"AS_T1_deplete_0": {'species': 'Bacteria', 'AS_type': 'Deplete human', 'Actual AS type': 'Deplete', 'Alt run name':'sample_0', 'library':7}}

def run():
    dfs=[]
    for run in runs:
        print(run)
        df=getFoldEnrichmentTime(f=f'/mnt/data/analysis/nick/adaptive_sampling/combined/results/{run}/mixdata.csv', 
                                species=runs[run]['species'], AS_type=runs[run]['AS_type'])
        df['run']=run
        df['species']=runs[run]['species']
        df['AS type']=runs[run]['Actual AS type']
        df['Alt run name']=runs[run]['Alt run name']
        df['library']=runs[run]['library']
        dfs.append(df)

    df=pd.concat(dfs)
    df['time']=pd.to_datetime(df['start_time_AS'])
    df['start_time']=df.groupby('run')['time'].transform(min)
    AS_S1_Ecoli_10_4_enrich_0_END_TIME=df[df['run']=='AS_S1_Ecoli_10_4_enrich_0']['time'].max()
    time_delta = df['start_time'] - AS_S1_Ecoli_10_4_enrich_0_END_TIME

    df['run time']=np.where(df['run']=='AS_S1_Ecoli_10_4_deplete_0', df['time']-(df['start_time']+time_delta), df['time']-df['start_time'])
    df['minutes']=df['run time']/60 
    df['hours']=df['minutes']/60
    df.to_csv('fold_enrichment.csv',index=False)
    return df

def plot(df):
    #df['time']=pd.to_datetime(df['start_time_AS'])
    #df['start_time']=df.groupby('run')['time'].transform(min)
    #AS_S1_Ecoli_10_4_enrich_0_END_TIME=df[df['run']=='AS_S1_Ecoli_10_4_enrich_0']['time'].max()
    #time_delta = df['start_time'] - AS_S1_Ecoli_10_4_enrich_0_END_TIME

    #df['run time']=np.where(df['run']=='AS_S1_Ecoli_10_4_deplete_0', df['time']-(df['start_time']+time_delta), df['time']-df['start_time'])
    #df['minutes']=df['run time']/60 
    #df['hours']=df['minutes']/60
    df['hours']=pd.to_timedelta(df['hours'])
    
    ax=sns.lineplot(x='hours', y='Fold enrichment', data=df, hue='Alt run name', style='AS type')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    #plt.show()
    plt.savefig('fold_enrichment_over_time.pdf')

if __name__ == '__main__':
    df=run()
    df=pd.read_csv('fold_enrichment.csv')
    df['run time']=pd.to_timedelta(df['run time'])
    df['seconds']=df['run time'].dt.total_seconds()
    df.drop_duplicates(subset=['run','seconds'], inplace=True, keep='last')
    df.to_csv('fold_enrichment2.csv',index=False)
    plot(df)
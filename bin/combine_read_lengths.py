#!/usr/bin/env python3
import pandas as pd
import numpy as np

runs=['AS_S1_Ecoli_10_4_deplete_0',
'AS_S1_Ecoli_10_4_enrich_0',
'AS_S2_1_Ecol_10_4_deplete_0',
'AS_S2_1_Ecol_10_4_enrich_0',
'AS_S2_1_Saur_10_4_enrich_0',
'AS_S2_2_Ecol_10_4_deplete_0',
'AS_S2_2_Ecol_10_4_enrich_0',
'AS_S2_2_Saur_10_4_enrich_0',
'AS_Staph_50_50_deplete_0',
'AS_T1_deplete_0',
'saur50_human50_20k_ligation']

runs=['sample_1a', 'sample_2a', 'sample_3a', 'sample_4a', 'sample_6a', 'sample_7a', 'sample_PBS_neg1', 'AS_T1_deplete_0']

dfs=[]
for run in runs:
    df=pd.read_csv(f'/mnt/data/analysis/nick/adaptive_sampling/combined/results/{run}/read_lengths.csv')
    df['run']=run
    df['Species/human']=np.where(df['Read type'].isin(['Human','Other']),df['Read type'],'Species')
    df['AS/CONTROL']=np.where(df['Condition'] == 'control condtion 1','control','AS')
    df['Exp_type']=df['Species/human']+'_'+df['AS/CONTROL']
    dfs.append(df)

df=pd.concat(dfs)

# create multi index with run and exp_type
#df.set_index(['run','Exp_type'],inplace=True)

# pivot table
df=df[['run','Exp_type','count','mean','std','min','25%','50%','75%','max']]
df2=df.pivot(index="run", columns="Exp_type")
print(df2)
df2.to_csv('adaptive_sampling_read_lengths.csv')


dfs=[]
for run in runs:
    df=pd.read_csv(f'/mnt/data/analysis/nick/adaptive_sampling/combined/results/{run}/Yields.csv')
    df['run']=run
    df['Species/human']=np.where(df['Read type'].isin(['Human','Other']),df['Read type'],'Species')
    df['AS/CONTROL']=np.where(df['Group'] == 'control condtion 1','control','AS')
    df['Exp_type']=df['Species/human']+'_'+df['AS/CONTROL']
    df=df[['run','Exp_type','count','max']]
    dfs.append(df)

df=pd.concat(dfs)
df.rename(columns={'count':'num_reads','max':'bases'},inplace=True)
df2=df.pivot(index="run", columns="Exp_type")
print(df2)
df2.to_csv('adaptive_sampling_yields.csv')

dfs=[]
for run in runs:
   df=pd.read_csv(f'/mnt/data/analysis/nick/adaptive_sampling/combined/results/{run}/run_summary.tsv', sep='\t')
   df['run']=run
   df['Species/human']=np.where(df['Read type'].isin(['Human','Other']),df['Read type'],'Species')
   df['AS/CONTROL']=np.where(df['Group'] == 'control condtion 1','control','AS')
   df['Exp_type']=df['Species/human']+'_'+df['AS/CONTROL']
   dfs.append(df)
   df=df[['run','Exp_type','Condition total active channels']]

df=pd.concat(dfs)
df.rename(columns={'count':'num_reads','max':'bases'},inplace=True)
df2=df.pivot(index="run", columns="Exp_type")
print(df2)
df2.to_csv('adaptive_sampling_active_pores.csv')

dfs=[]
cols=['read_id','Read type', 'Condition', 'sequence_length_template', 'unblocked']
for run in runs:
    print(run)
    df=pd.read_csv(f'/mnt/data/analysis/nick/adaptive_sampling/combined/results/{run}/mixdata.csv', sep='\t')
    df=df[cols]
    df['run']=run
    df=df.groupby(['Read type', 'Condition']).sample(1000, replace=True)
    df.drop_duplicates(subset=['read_id'], inplace=True)
    df['Species/human']=np.where(df['Read type'].isin(['Human','Other']),df['Read type'],'Species')
    df['AS/CONTROL']=np.where(df['Condition'] == 'control condtion 1','control','AS')
    df['Exp_type']=df['Species/human']+'_'+df['AS/CONTROL']
    dfs.append(df)

df=pd.concat(dfs)
df.to_csv('adaptive_sampling_read_lengths_raw.csv')


dfs=[]
cols=['read_id','hours', 'Condition', 'Condition active channels', 'run time']
for run in runs:
    print(run)
    df=pd.read_csv(f'/mnt/data/analysis/nick/adaptive_sampling/combined/results/{run}/mixdata.csv', sep='\t')
    df=df[cols]
    df['run']=run
    df=df.groupby(['hours', 'Condition']).sample(60, replace=True)
    df['run time']=pd.to_timedelta(df['run time'])
    df['seconds']=df['run time'].dt.total_seconds()
    df['AS_CONTROL']=np.where(df['Condition'] == 'control condtion 1','control','AS')
    df.drop_duplicates(subset=['read_id'], inplace=True)
    dfs.append(df)

df=pd.concat(dfs)
df.to_csv('adaptive_sampling_active_pores_sampled.csv')
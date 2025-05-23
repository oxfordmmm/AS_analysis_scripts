#!/usr/bin/env python3
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")
import sys
import gzip
from argparse import ArgumentParser, SUPPRESS
from ete3 import NCBITaxa
ncbi = NCBITaxa()

df=pd.read_csv(sys.argv[1], sep='\t')
df['Is Human read']=np.where(df['taxID']==9606, True, False)
#df['Group']=np.where(df['Condition']=='Condtion0', 'AS human deplete', 'Control')
df['Group']=df['Condition']

# is species read?
species=str(sys.argv[4])
taxID=int(sys.argv[5])
species_taxids = ncbi.get_descendant_taxa(taxID, intermediate_nodes=True)
species_taxids.append(taxID)
df['Is {0} read'.format(species)]=np.where(df['taxID'].isin(species_taxids), True, False)

# read type
def read_type(r):
	if r['Is {0} read'.format(species)] ==True:
		r['Read type']= species
	elif r['Is Human read'] == True:
		r['Read type']= 'Human'
	else:
		r['Read type']= 'Other'
	return r

df['Read type']=np.where(df['taxID'].isin(species_taxids), species, 'Other')
df['Read type']=np.where(df['taxID']==9606, 'Human', df['Read type'])
#df['Read type']=np.where(df['Read type']==False, df['Read type'],'Other')
#df=df.apply(read_type, axis=1)


# culmulative bases
df['time']=pd.to_datetime(df['start_time'])
start_time=df['time'].min()
df['run time']=df['time']-start_time
df['minutes']=df['run time']/60 
df['hours']=df['minutes']/60
df.sort_values(by='run time',inplace=True)
df['Group cumulative bases']=df.groupby(['Group','Read type'])['sequence_length_template'].transform('cumsum')

# active channels
df['10mins']=df['minutes']=df['minutes']/10
df['mux channel']=df['mux'].map(str) + df['channel'].map(str)
df['Condition active channels']=df.groupby(['Condition','10mins'])['mux channel'].transform('nunique')
df['Condition total active channels']=df.groupby(['Condition'])['mux channel'].transform('nunique')
df['Condition active channel time']=df.groupby(['Group'])['duration'].transform('cumsum')

# bases per active channels BPAC
df['BPAC']=df['Group cumulative bases']/df['Condition total active channels']

# prop reads
df['5mins']=df['minutes']=df['minutes']/5
df['Group reads per 5mins']=df.groupby(['Group','Read type','5mins'])[['read_id']].transform('count')
df['Condtion reads per 5mins']=df.groupby(['Group','5mins'])[['read_id']].transform('count')
df['Group % of reads sequenced']=(df['Group reads per 5mins']/df['Condtion reads per 5mins'])*100

## read length boxplots
def read_length_boxplot(df):
	g1=sns.boxplot(x='Read type',
		y='sequence_length_template',
		hue='Group',
		data=df)

	#g1.set_yscale('log')
	g1.set(ylim=(0, 15000))
	plt.savefig('read_length_boxplot.pdf')
	plt.clf()

	g=df.groupby(['Read type', 'Condition'])['sequence_length_template'].describe()
	print(g)
	g.to_csv('read_lengths.csv')

def read_length_violinplot(df):
	g1=sns.violinplot(x='Read type',
		y='sequence_length_template',
		hue='Group',
		data=df)

	#g1.set_yscale('log')
	g1.set(ylim=(0, 15000))
	plt.savefig('read_length_violinplot.pdf')
	plt.clf()


## yield over time
def yield_over_time(df):
	
	g2=sns.lineplot(x='hours', y='Group cumulative bases', 
			data=df,
			hue='Group',
			style='Read type')
	#g2.set(yscale='log')
	for l in g2.lines:
		y = l.get_ydata()
		#print(y)
		if len(y)>0:
			v=y[-1]/1000000
			print(f'{v:.2f}MB')
			g2.annotate(f'{v:.2f}MB', xy=(1,y[-1]), xycoords=('axes fraction', 'data'), 
						ha='left', va='center', color=l.get_color())

	plt.tight_layout()
	plt.savefig('yield_time.pdf')
	plt.clf()
	return df

def yield_over_time_panel(df):
	g21=sns.relplot(x='hours', y='Group cumulative bases', 
			data=df,
			hue='Group',
			style='Read type',
			col='Read type',
			row='Group',
			kind="line",
			legend = False)

	#g2.set(yscale='log')
	for ax in g21.axes.ravel():
		for l in ax.lines:
		    y = l.get_ydata()
		    #print(y)
		    if len(y)>0:
			    v=y[-1]/1000000
			    print(f'{v:.2f}MB')
			    ax.annotate(f'{v:.2f}MB', xy=(1,y[-1]), xycoords=('axes fraction', 'data'), 
						ha='left', va='center', color=l.get_color())

	plt.tight_layout()
	plt.savefig('yield_time_panel.pdf')
	plt.clf()
	return df

def bpac_over_time(df):
	
	g21=sns.lineplot(x='hours', y='BPAC', 
			data=df,
			hue='Group',
			style='Read type')
	#g2.set(yscale='log')
	for l in g21.lines:
		y = l.get_ydata()
		#print(y)
		if len(y)>0:
			v=y[-1]/1000000
			print(f'{v:.2f}MB')
			g21.annotate(f'{v:.2f}MB', xy=(1,y[-1]), xycoords=('axes fraction', 'data'), 
						ha='left', va='center', color=l.get_color())

	plt.tight_layout()
	plt.savefig('bpac_time.pdf')
	plt.clf()
	return df

def bpac_over_time_panel(df):
	g21=sns.relplot(x='hours', y='BPAC', 
			data=df,
			hue='Group',
			style='Read type',
			col='Read type',
			row='Group',
			kind="line",
			legend = False)

	#g2.set(yscale='log')
	for ax in g21.axes.ravel():
		for l in ax.lines:
		    y = l.get_ydata()
		    #print(y)
		    if len(y)>0:
			    v=y[-1]/1000000
			    print(f'{v:.2f}MB')
			    ax.annotate(f'{v:.2f}MB', xy=(1,y[-1]), xycoords=('axes fraction', 'data'), 
						ha='left', va='center', color=l.get_color())

	plt.tight_layout()
	plt.savefig('bpac_time_panel.pdf')
	plt.clf()
	return df

## ## read length over time
def read_length_over_time(df):
	g3=sns.relplot(x='hours',y='sequence_length_template',
			data=df,
			hue='Group',
			style='Read type',
			col='Read type',
			row='Group',
			kind="line",
			legend = False)
	plt.tight_layout()
	plt.savefig('read_length_over_time.pdf')
	plt.clf()


## number or prop or reads over time
def prop_reads_over_time(df):
	# plot
	g4=sns.lineplot(x='hours',
				y='Group % of reads sequenced',
				data=df,
				hue='Group',
				style='Read type')

	plt.tight_layout()
	plt.savefig('prop_reads_over_time.pdf')
	plt.clf()

def prop_reads_over_time_panel(df):
	# plot
	g4=sns.relplot(x='hours',
				y='Group % of reads sequenced',
				data=df,
				hue='Group',
				style='Read type',
				col='Read type',
				row='Group',
				kind="line",
				legend = False)

	plt.tight_layout()
	plt.savefig('prop_reads_over_time_panel.pdf')
	plt.clf()

## prop channels over time
def active_chanels(df):
	#df['10mins']=df['minutes']=df['minutes']/10
	#df['Condition active channels']=df.groupby(['Condition','10mins'])['channel'].transform('nunique')

	g5=sns.lineplot(x='hours',
					y='Condition active channels',
					data=df,
					hue='Group')
	
	plt.tight_layout()
	plt.savefig('channels_over_time.pdf')
	plt.clf()

def active_chanels_time(df):
	#df['10mins']=df['minutes']=df['minutes']/10
	#df['Condition active channels']=df.groupby(['Condition','10mins'])['channel'].transform('nunique')

	g5=sns.lineplot(x='hours',
					y='Condition active channel time',
					data=df,
					hue='Group')
	
	plt.tight_layout()
	plt.savefig('channels_time_over_time.pdf')
	plt.clf()

read_length_boxplot(df)
read_length_violinplot(df)
df=yield_over_time(df)
read_length_over_time(df)
prop_reads_over_time(df)
active_chanels(df)
bpac_over_time(df)
bpac_over_time_panel(df)
prop_reads_over_time_panel(df)
active_chanels_time(df)
yield_over_time_panel(df)

## yields
g=df.groupby(['Read type','Group'])['Group cumulative bases'].describe()
print(g)
g.to_csv('Yields.csv')

# whole run summary
g1=df.groupby(['Read type','Group'])['Group cumulative bases'].max().reset_index() 
g2=df.groupby(['Read type','Group'])['Condition total active channels'].max().reset_index() 
g3=df.groupby(['Read type', 'Group'])[['read_id']].count().reset_index() 
g4=g1.merge(g2, on=['Read type','Group'], how='left') 
g4=g4.merge(g3, on=['Read type','Group'], how='left') 
g4['Pan condition bases']=g4.groupby(['Group'])[['Group cumulative bases']].transform(sum)
g4['Group bases yield %']=(g4['Group cumulative bases']/g4['Pan condition bases'])*100
g4['Mean read length']=df.groupby(['Read type', 'Condition'])['sequence_length_template'].transform('mean')
g4['Condition diff']=g4.groupby(['Group'])['Group bases yield %'].transform('diff')
g4['Read type flowcell bases']=g4.groupby(['Read type'])[['Group cumulative bases']].transform(sum)


g4['Total group reads']=g4.groupby(['Group'])['read_id'].transform(sum)
g4['Read % in group']=100*(g4['read_id']/g4['Total group reads'])

g4['Run name']=sys.argv[3]
g4.to_csv('run_summary.tsv',index=False,sep='\t')

## unblocked reads
unblocked=open(sys.argv[2],'rt').read().split('\n')
df['unblocked']=np.where(df['read_id'].isin(unblocked), True, False)

g=df.groupby(['Read type', 'Group','unblocked'])[['read_id']].count()
g.reset_index(inplace=True)
print(g)
g['group total']=g.groupby(['Read type', 'Group'])[['read_id']].transform(sum)
g['group unblock %']=100*(g['read_id']/g['group total'])

g.rename(columns={'read_id':'read numbers'}, inplace=True)
print(g)

g.to_csv('Group_counts.csv')

df.to_csv('mixdata.csv',sep='\t',index=False)


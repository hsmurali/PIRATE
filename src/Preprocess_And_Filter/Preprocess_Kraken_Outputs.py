#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from os import listdir

metacarvel_path = '/Users/harihara/Mount-3/hmp_scaffolds/stool/'
kraken_path = '/Users/harihara/Mount-3/hmp_scaffolds/stool/kraken_output/'
virsorter_path = '/Users/harihara/Research-Activities/Data/Virsorter/Stool/'

motif_types = ['four_bubbles','three_bubbles','complex_bubbles','repeats','tandem_repeats']

samples = listdir(metacarvel_path)
samples.sort()

def Get_Contig_List(filepath):
    lines = open(filepath).readlines()
    op_list = []
    for line in lines:
        line = line.replace('\n',"")
        line = list(set(line.split('\t')))
        op_list += line
    return op_list

def Parse_Phylogenetic_Classification(phyla):
    phyla = phyla.split(';')
    splits = [] 
    for x in phyla: 
        if ("/" not in x and "group" not in x and "sub" not in x and "Sub" not in x):
            splits.append(x)
    return len(splits)

def Parse_Virsorter_Outputs(filepath):
    df = pd.read_csv(filepath+'viral_contigs_in_features.txt', delimiter = '\t')
    df['CONTIG'] = df['CONTIG'].astype('str')
    df = df.rename(columns = {'CONTIG':'ContigID','SAMPLE':'SampleID'})
    return df[['SampleID','ContigID','CATEGORY', 'FEATURE']]




# In[ ]:


Outputs = dict()
counter = 0

for sample in samples:
    if sample.startswith('SRS'):
        counter += 1
        motifs_path = metacarvel_path + sample + '/' + sample + '_scaffolds/'
        sample_kraken_path = kraken_path+'/'+sample+'_kraken.labels'
        df_kraken = pd.read_csv(sample_kraken_path, names = ['ContigID','Phylum-Label'], sep = '\t')
        df_kraken.set_index('ContigID', inplace = True)
        print (sample, len(df_kraken))
        for motif in motif_types:
            contigs = Get_Contig_List(motifs_path+motif)
            try:
                df = pd.DataFrame(data = {'ContigID':contigs})
                df.loc[:,'SampleID'] = sample
                df.set_index('ContigID', inplace = True)
                df = df_kraken.join(df, how = 'right')
                df = df.reset_index()
                df['Phylum-Label'] = df['Phylum-Label'].fillna('Unclassified')
                df['Taxon-Label'] = df['Phylum-Label'].apply(Parse_Phylogenetic_Classification)
            except ValueError:
                df = pd.DataFrame()
            try:
                Outputs[motif] = Outputs[motif].append(df, ignore_index = True)
            except KeyError:
                Outputs[motif] = df       

df_stool = pd.DataFrame()
for motif in motif_types:
    df_stool = df_stool.append(Outputs[motif][['SampleID','ContigID','Phylum-Label','Taxon-Label']])


# In[ ]:


op_path = '/Users/harihara/Research-Activities/Data/Kraken_Updated/'
df_stool.to_pickle(op_path+'Stool.pkl')


# In[ ]:


df_virsorter = Parse_Virsorter_Outputs(virsorter_path)
df_virsorter.to_pickle('/Users/harihara/Research-Activities/Data/Virsorter/Stool/Stool.pkl')


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from functools import partial
from os import listdir
global counter

counter = 0
data_dir = '/Users/harihara/Research-Activities/Data/'


# In[ ]:


def Aggregate_Fraction_of_Coverage(group):
    global counter
    counter += 1
    print(counter)
    contig_length = group.iloc[0]['Length']
    idx = (group.groupby(['Q_Align_Start'])['Perfect_Alignment_Length'].transform(max) == 
           group['Perfect_Alignment_Length'])
    group = group[idx]
    group = group.reset_index()
    group = group.sort_values(by = ['Q_Align_Start','Q_Align_End'])
   
    gaps = np.array(group['Num_Gap_Openings'])
    start = np.array(group['Q_Align_Start'])
    end = np.array(group['Q_Align_End'])
    align_lengths = np.array(group['Alignment_Length'])
    
    fraction_covered = align_lengths[0] - gaps[0]
    start_1 = start[0]
    end_1 = end[0]
    
    for i in range(1, len(group)):
        start_2 = start[i]
        end_2 = end[i]
        alignment_length = (align_lengths[i]) - (gaps[i])
        if start_2 >= end_1:
            start_1 = start_2
            end_1 = end_2
            fraction_covered += (alignment_length)
        elif start_2 >= start_1 and start_2 < end_1:
            if end_2 > end_1:
                alignment_length = end_2 - end_1 - gaps[i] - gaps[i-1]
                fraction_covered += alignment_length
                end_1 = end_2
    return fraction_covered*100.0/contig_length

def Get_Blast_Stats(filepath): 
    global counter
    counter = 0
    df_pkl = pd.read_pickle(filepath)
    df_pkl['Perfect_Alignment_Length'] = (df_pkl['Percent-Identity']*df_pkl['Alignment_Length'] - 
                                          df_pkl['Num_Gap_Openings'])
    df_pkl = df_pkl.sort_values(by = ['Q_Align_Start','Q_Align_End'])
    df_pkl['Hits'] = 0
    df_grouped = df_pkl.groupby(['SampleID','ContigID','Length'])
    print(len(df_grouped))
    df_alignment_stats = df_grouped['Perfect_Alignment_Length'].agg(['max','mean'])
    df_alignment_stats.rename(columns = {'max':'Max_Perfect_Alignment','mean':'Mean_Perfect_Alignment'}, inplace = True)
    df_E_Val_Stats = df_grouped['E-Val'].agg(['max','mean'])
    df_E_Val_Stats.rename(columns = {'max':'Max_E_Val','mean':'Mean_E_Val'}, inplace = True)
    df_Hits = df_grouped['Hits'].count()
    df_stats = pd.DataFrame()
    df_stats = df_alignment_stats
    df_stats = df_stats.join(df_E_Val_Stats)
    df_stats = df_stats.join(df_Hits)
    df_stats['Breadth_Coverage'] = df_grouped.apply(Aggregate_Fraction_of_Coverage)
    return df_stats


# In[ ]:


in_filepath = data_dir+'Blast/Stool/PKL/'
op_stats = data_dir+'Blast/Stool/Blast_Stats/'
files = listdir(in_filepath)
for f in files:
    print(f)
    if f[0] != '.':
        df_stats = Get_Blast_Stats(in_filepath+f)
        df_stats.to_pickle(op_stats+f)


# In[ ]:





# In[ ]:





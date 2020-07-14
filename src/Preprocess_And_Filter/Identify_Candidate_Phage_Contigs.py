#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import listdir


# In[ ]:


rcParams = {'font.size': 25 , 'font.weight': 'normal', 'font.family': 'sans-serif',
            'axes.unicode_minus':False, 'axes.labelweight':'normal'}
plt.rcParams.update(rcParams)


# In[ ]:


data_dir = '/Users/harihara/Research-Activities/Data/'
bubble_filepath = data_dir+'Bubble_Length_Statistics_Updated_Nomenclature/stool/'
repeats_path = data_dir+'Repeat_Statistics/stool/'
blast_path = data_dir+'Blast/Stool/Blast_Stats/'
Kraken_Path = data_dir+'Kraken_Updated/Stool.pkl'
virsorter_path = data_dir+'Virsorter/Stool/Stool.pkl'


# In[ ]:


def Return_Blast_Outputs(filepath):
    files = listdir(filepath)
    df_blast = pd.DataFrame()
    for f in files:
        df = pd.read_pickle(filepath+f)
        df = df.drop_duplicates()
        df_blast = df_blast.append(df, ignore_index = False)
    df_blast = df_blast.drop_duplicates().reset_index('Length')
    return df_blast

def Fill_Missing_Lengths(filepath, df, min_contig_length = 5000):
    df_lengths = pd.read_pickle(filepath)
    del df['Length']
    idx = df.index
    df_lengths_select = df_lengths.loc[idx,:]
    df = df.join(df_lengths_select, how = 'left')
    df = df.drop_duplicates()
    df = df[(df['Motif_Type'] != 'Interspersed') | (df['Length'] > min_contig_length)]
    return df

def Return_All_Feature_Contigs(path_length=0):
    files = listdir(bubble_filepath)
    files.sort()
    contigs_data, repeats_data = [], []
     
    for f in files:
        print(f)
        if f[0] != 'S':
            continue
        sample_id = f[:-4]
        df_repeats = pd.read_csv(repeats_path+f)
        df_repeats = df_repeats.rename(columns = {'Repeat_Type':'Motif_Type','Contig':'ContigID'})
        df_repeats.loc[:,'SampleID'] = sample_id
        df_repeats = df_repeats[['SampleID','ContigID','Motif_Type']]
        repeats_data += (list(df_repeats.T.to_dict().values()))
        df = pd.read_csv(bubble_filepath+f, low_memory = False)
        length_cols = [col for col in df.columns if 'Mean_' in col]
        df_length = df[length_cols]
        df_length_gr_11000 = df_length > path_length
        new_col = df_length_gr_11000.sum(axis = 1)
        df.loc[:,'Count_Gr_11000'] = new_col
        df = df[df['Count_Gr_11000'] > 0]
        for i in range(0, len(df)):
            row = df.iloc[i].to_dict()
            for m in length_cols:
                if row[m] > path_length:
                    p = m.index('_')
                    contigs = row['Path_'+m[p+1:]].split(' ')
                    cols = ['SampleID']*len(contigs) + ['ContigID']*len(contigs) + ['Motif_Type']*len(contigs)
                    values = [sample_id]*len(contigs) + contigs + [row['Bubble_Type']]*len(contigs)
                    contigs_data.append(dict(zip(cols, values)))
    df_contigs = pd.DataFrame(contigs_data)
    df_contigs = df_contigs.drop_duplicates()
    df_repeats = pd.DataFrame(repeats_data)
    df_contigs = df_contigs.append(df_repeats, ignore_index = True, sort = True)
    df_contigs = df_contigs.set_index(['SampleID','ContigID'])
    return df_contigs


# In[ ]:


df_blast = Return_Blast_Outputs(blast_path)
df_kraken = pd.read_pickle(Kraken_Path).set_index( ['SampleID','ContigID'])[['Taxon-Label']]
df_virsorter = pd.read_pickle(virsorter_path).set_index(
    ['SampleID','ContigID']).rename(
    columns = {'CATEGORY':'Virsorter-Category'})


# In[ ]:


df_contigs_potential_phages = Return_All_Feature_Contigs()
df_contigs_potential_phages['Contigs_of_Interest'] = True
df_all_results = df_blast.join(df_kraken, how = 'outer').join(
    df_virsorter, how = 'outer').join(df_contigs_potential_phages, how = 'right')
df_all_results[['Breadth_Coverage']] = df_all_results[['Breadth_Coverage']].fillna(0)
df_all_results[['Virsorter-Category']] = df_all_results[['Virsorter-Category']].fillna("NA")
df_all_results[['Taxon-Label']] = df_all_results[['Taxon-Label']].fillna(-1)
df_all_results = Fill_Missing_Lengths(data_dir+'Stool_Contig_Length.pkl', 
                                      df_all_results)
del df_all_results['FEATURE'], df_all_results['Unnamed: 0']
df_all_results.loc[df_all_results['Taxon-Label'] >= 10, 'Taxon-Label'] = 10
df_potential_phages = df_all_results[df_all_results['Breadth_Coverage'] < 50]
df_potential_phages = df_potential_phages[df_potential_phages['Taxon-Label'] < 7]

print (p, len(df_potential_phages))


# In[ ]:


fig_hist, ax_hist = plt.subplots(1,1,figsize = (16,10))
df_all_results['Breadth_Coverage'].hist(density = True, figsize = (16,10), 
                                        cumulative = True, histtype = 'step',
                                        linewidth = 4, ax = ax_hist, bins = 50)
ax_hist.set_xlabel('Breadth of Coverage')
fig_hist.tight_layout()


# In[ ]:


fig_bar, ax_bar = plt.subplots(1,1,figsize = (16,10))
df_all_results['Count'] = 0
df_all_grouped = df_all_results.groupby('Taxon-Label').count()['Count']
df_all_grouped.plot.bar(legend = False, ax = ax_bar)
ax_bar.set_xlabel('Taxanomic Label')
ax_bar.set_ylabel('No.of Contigs')
ax_bar.set_xticklabels(['Unclassified/Root','Domain','Kingdom','Phylum','Class',
                        'Order','Family','Genus','Species','Strains'],rotation=90)
fig_bar.tight_layout()


# In[ ]:


fig_vir_bar, ax_vir_bar = plt.subplots(1,1,figsize = (16,10))
df_all_results['Count'] = 0
df_all_grouped = df_all_results.groupby('Virsorter-Category').count()['Count']
df_all_grouped.plot.bar(legend = False, ax = ax_vir_bar)
ax_vir_bar.set_xlabel('Virsorter-Category')
ax_vir_bar.set_ylabel('No.of Contigs')
fig_vir_bar.tight_layout()


# In[ ]:


df_potential_phages.to_csv(data_dir+ 'Potential_Phages/Filtered_Phage_List.csv')


# In[ ]:





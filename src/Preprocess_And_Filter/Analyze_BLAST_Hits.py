#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from os import listdir
import matplotlib.pyplot as plt


# In[ ]:


rcParams = {'font.size': 25 , 'font.weight': 'normal', 'font.family': 'sans-serif',
            'axes.unicode_minus':False, 'axes.labelweight':'normal'}
plt.rcParams.update(rcParams)


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


# In[ ]:


df_blast = Return_Blast_Outputs('/Users/harihara/Research-Activities/Data/Blast/Stool/Blast_Stats/')
plots_path = '/Users/harihara/Research-Activities/Plots/Potential_Phages/'


# In[ ]:


fig_hist, ax_hist = plt.subplots(1,1,figsize = (16,10))
df_blast[['Breadth_Coverage']].hist(bins = 50, density = True, cumulative = True, color = 'blue',
                                    histtype = 'step', linewidth = 3, ax = ax_hist)
ax_hist.set_xlabel('Breadth of Coverage (in %)')
ax_hist.set_ylabel('Fraction of Contigs')
ax_hist.set_title('')
fig_hist.tight_layout()
fig_hist.savefig(plots_path+'CDF-Breadth-Coverage-Blast.pdf')


# In[ ]:


fig_e, ax_e = plt.subplots(1,2,figsize=(16,10))
ax_e[0].boxplot(df_blast['Mean_E_Val'], showfliers = False,
                boxprops = {'color':'blue', 'linewidth' : 4}, 
                medianprops = {'color' : 'red', 'linewidth':4},
                whiskerprops = {'color' : 'blue'})
ax_e[0].set_xticklabels('')
ax_e[0].set_xlabel('Average E-Value')

df_blast.plot.scatter('Breadth_Coverage','Mean_E_Val', color = 'blue', s = 50, ax = ax_e[1])
ax_e[1].set_ylabel('Average E-Value')
ax_e[1].set_xlabel('Breadth of Coverage')
ax_e[1].set_ylim([0,3])
fig_e.tight_layout()

fig_e.savefig(plots_path+'E-Val-Stats-Blast.pdf')


# In[ ]:


fig_length, ax_length = plt.subplots(1,1,figsize = (16,10))
df_blast.plot.scatter('Breadth_Coverage','Length', color = 'blue', s = 50, ax = ax_length)
ax_length.set_xlabel('Breadth of Coverage')
ax_length.set_ylabel('Contig Length')
fig_length.tight_layout()
fig_length.savefig(plots_path+'Length-Breadth-Coverage-Blast.pdf')


# In[ ]:





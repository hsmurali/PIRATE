#!/usr/bin/env python
# coding: utf-8

# In[48]:


import numpy as np
import pandas as pd

data_path = '/Users/harihara/Research-Activities/Data/'

def Prepare_Metadata_Uniprot_Data(filepath):
    filebuf = open(filepath, "r").readlines()
    op = []
    for line in filebuf:
        if line[0] == '>':
            splits = line.split(' ')
            d = dict()
            d['Phage_ID'] = splits[0][1:]
            ind_start = 1
            for i in range(1,len(splits)):
                if 'OS=' in splits[i]:
                    d['Protein Identifier'] = " ".join(splits[1:i])
                    ind_start = i
                if 'OX=' in splits[i]:
                    ind_end = i
                    d['Organism Name'] = " ".join(splits[ind_start:ind_end])
            op.append(d)
    df = pd.DataFrame(op)
    return df

def Integrase():
    filepath = data_path + 'uniprot-organism_phage+integrase.fasta'
    op = []
    for f in open(filepath,'r').readlines():
        if f[0] == '>':
            d = {'Phage_ID':f.split(' ')[0][1:], 'Integrase':1} 
            op.append(d)
    df_integrase = pd.DataFrame(op)
    return df_integrase


# In[49]:


df_Meta_Uniprot = Prepare_Metadata_Uniprot_Data(data_path+'uniprot-organism_phage.fasta')
df_Integrase = Integrase()
df_Meta_Uniprot.set_index('Phage_ID', inplace = True)
df_Integrase.set_index('Phage_ID', inplace = True)
df_Meta_Uniprot = df_Meta_Uniprot.join(df_Integrase, how = 'left')
df_Meta_Uniprot['Integrase'].fillna(0, inplace = True)


# In[50]:


df_Meta_Uniprot.to_csv(data_path+'Uniprot_Meta_Data.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





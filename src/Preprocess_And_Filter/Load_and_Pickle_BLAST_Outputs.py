#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from os.path import isdir
from os import mkdir
from functools import partial


# In[ ]:


def Extract_sample_contigs(Query, options, repeat):
    split = Query.split('_')
    if options == 'S':
        return split[0]
    if options == 'C':
        if repeat == True:
            return "_".join(split[2:])
        if repeat == False:
            return "_".join(split[3:])

def Process_Chunk(chunk, repeat):
    chunk['SampleID'] = chunk['Query'].apply(partial(Extract_sample_contigs, 
                                               options = 'S', repeat = repeat))
    chunk['ContigID'] = chunk['Query'].apply(partial(Extract_sample_contigs, 
                                           options = 'C', repeat = repeat))
    chunk['ContigID'] = chunk['ContigID'].astype('str')
    chunk = chunk.set_index(['SampleID','ContigID'])
    return chunk

def Make_Pickle_Blast_Outputs(filepath, repeat):
    chunks = pd.read_csv(filepath, delimiter = '\t', 
                         names =  ['Query','Subject','Percent-Identity',
                                   'Alignment_Length','Num_Mismatch', 
                                   'Num_Gap_Openings','Q_Align_Start',
                                   'Q_Align_End', 'S_Align_Start', 
                                   'S_Align_End', 'E-Val','Bit-Score'], 
                         dtype = dict(zip(['Query','Subject','Percent-Identity',
                                           'Alignment_Length','Num_Mismatch', 
                                           'Num_Gap_Openings','Q_Align_Start',
                                           'Q_Align_End', 'S_Align_Start', 
                                           'S_Align_End', 'E-Val','Bit-Score'],
                                          ['str','str','float','int','int','int',
                                           'int','int','int','float','float'])),
                         low_memory = False, chunksize = 500000)
    df = pd.DataFrame()
    ctr = 0
    for chunk in chunks:
        print(ctr)
        ctr+=1
        chunk = Process_Chunk(chunk, repeat)
        df = df.append(chunk, ignore_index = False)
    return df

def Load_Contig_Lengths_File(filepath):
    df_lengths = pd.read_csv(filepath, index_col = ['SampleID','ContigID'])
    df_lengths.to_pickle(filepath[:-4]+".pkl")
    return df_lengths


# In[ ]:


contig_length_filepath = '/Users/harihara/Research-Activities/Data/Stool_Contig_Length.csv'
blast_path = '/Users/harihara/Research-Activities/Data/Blast/Stool/stool_'
op_path = '/Users/harihara/Research-Activities/Data/Blast/Stool/PKL/'
if not isdir(op_path):
    mkdir(op_path)


# In[ ]:


df_length = Load_Contig_Lengths_File(contig_length_filepath)
motif_types = ['three_bubble_contigs','four_bubble_contigs','complex_bubble_contigs',
               'interspersed_repeats','tandem_repeats']
for motif_type in motif_types[:3]:
    print(motif_type)
    if 'repeats' in motif_type:
        repeats = True
    else:
        repeats = False
    blast_motif_path = blast_path+motif_type+'.blast.out'
    print(blast_motif_path)
    df_blast = Make_Pickle_Blast_Outputs(blast_motif_path, repeats)
    df_blast = df_blast.join(df_length, how = 'left')
    df_blast.to_pickle(op_path+'stool_'+motif_type+'.pkl')


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[ ]:


rcParams = {'font.size': 30 , 'font.weight': 'normal', 'font.family': 'sans-serif',
            'axes.unicode_minus':False, 'axes.labelweight':'normal'}
plt.rcParams.update(rcParams)


# In[ ]:


def Split_Contigs_ID(s, options):
    splits = s.split('_')
    if options == 'S':
        return splits[0]
    else:
        return "_".join(splits[1:])


# In[ ]:


df_lengths = pd.read_pickle('/Users/harihara/Research-Activities/Data/Stool_Contig_Length.pkl')


# In[ ]:


df_cluster = pd.read_csv('/Users/harihara/Mount/Potential_Phages_Updated/Clustering_Results.csv')
df_cluster['SampleID'] = df_cluster['Contigs'].apply(partial(Split_Contigs_ID, options = 'S'))
df_cluster['ContigID'] = df_cluster['Contigs'].apply(partial(Split_Contigs_ID, options = 'C'))
df_cluster = df_cluster.set_index(['SampleID','ContigID'])
indices = df_cluster.index.tolist()
df_lengths = df_lengths.loc[indices]
df_cluster = df_cluster.join(df_lengths)
del df_cluster['Contigs'], df_cluster['Unnamed: 0']


# In[ ]:


df_cluster_counts = df_cluster.max(axis = 0)

single = df_cluster_counts.loc[df_cluster_counts.index.str.startswith("Single_")]
single.index = single.index.str.replace("Single_","")
single.index = single.index.map(float)
single = single.sort_index()

complete = df_cluster_counts.loc[df_cluster_counts.index.str.startswith("Complete_")]
complete.index = complete.index.str.replace("Complete_","")
complete.index = complete.index.map(float)
complete = complete.sort_index()

average = df_cluster_counts.loc[df_cluster_counts.index.str.startswith("Average_")]
average.index = average.index.str.replace("Average_","")
average.index = average.index.map(float)
average = average.sort_index()

#ward = df_cluster_counts.loc[df_cluster_counts.index.str.startswith("Ward_")]
#ward.index = ward.index.str.replace("Ward_","")
#ward.index = ward.index.map(float)
#ward = ward.sort_index()

weighted = df_cluster_counts.loc[df_cluster_counts.index.str.startswith("Weighted_")]
weighted.index = weighted.index.str.replace("Weighted_","")
weighted.index = weighted.index.map(float)
weighted = weighted.sort_index()

centroid = df_cluster_counts.loc[df_cluster_counts.index.str.startswith("Centroid_")]
centroid.index = centroid.index.str.replace("Centroid_","")
centroid.index = centroid.index.map(float)
centroid = centroid.sort_index()

df_cluster_counts = pd.concat([single, complete, average], axis=1).rename(columns = 
                               {0:'Single', 1:'Complete', 2:'Average'}) 
                                #3:'Ward', 
#                                3:'Centroid',4:'Weighted'})

fig_cluster_counts, ax_cluster_counts = plt.subplots(1,1,figsize = (16,10))
df_cluster_counts.plot(marker = 'o', ms = 10, linewidth = 2.3, 
                       ax = ax_cluster_counts, linestyle = '-.')
ax_cluster_counts.scatter(0.61, df_cluster_counts.loc[0.61]['Complete'], marker = 'o', s = 500,
                          facecolors='none', edgecolors='indigo', linewidth = 2)
ax_cluster_counts.legend(ncol = 1)
ax_cluster_counts.set_xlabel('Cophenetic Distance')
ax_cluster_counts.set_ylabel('Number of Clusters')
ax_cluster_counts.set_xlim([0.35, 1.0])
ax_cluster_counts.grid()

fig_cluster_counts.tight_layout()
fig_cluster_counts.savefig('Cluster_Counts_Methods.pdf')


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (16,10))
ax.plot(complete.values, complete.index, linewidth = 5, color = 'blue')
ax.plot([0,complete[0.61]], [0.61, 0.61], color = 'black')
ax.plot([complete[0.61], complete[0.61]], [0.61, 0.4], color = 'black')
ax.scatter([complete[0.61]], [0.61], color = 'black', s = 150)
ax.set_xlim([0, 10000])
ax.set_ylim([0.4, 1])

ax.set_xlabel('Number of Clusters')
ax.set_ylabel('Cophenetic Distance')
fig.tight_layout()
fig.savefig('Cluster_Counts_Cophenetic_Distance.pdf')


# In[ ]:


columns_ward = df_cluster.columns[df_cluster.columns.str.startswith("Complete")].tolist()
Ward_Dict = {}
sample_diversity_dict = []
for c in columns_ward:
    Distance = c.split('_')[1]
    Method = c.split('_')[0]
    counts = df_cluster[[c]].reset_index().groupby(c).count()['SampleID'].tolist()
    Ward_Dict[float(Distance)] = counts


# In[ ]:


fig_counts, ax_counts = plt.subplots(1,1,figsize = (16,10))
keys = list(Ward_Dict.keys())
keys.sort()
Values = []
for k in keys:
    Values.append(Ward_Dict[k])
    
bp = ax_counts.boxplot(Values, whis = [5,90], showfliers = False, 
                  boxprops = {'color':'blue', 'linewidth':4},
                  medianprops = {'color':'orange', 'linewidth':3})
ax_counts.set_xticklabels(keys, rotation = 90)
ax_counts.set_xlabel('Cophenetic Distance Cutoff')
ax_counts.set_ylabel('# of Contigs in each cluster')

ax_counts_2 = ax_counts.twinx()
ctr = 4.5
ax_counts_2.plot([-2000]+complete.values.tolist(), marker = 'o', ms = 12.5, linewidth = 0.5, color = 'teal')
for x in complete.values.tolist()[4:28]:
    ax_counts_2.annotate(str(x), (ctr+0.35, x+2250), rotation = 50, fontsize = 27, color = 'red', alpha = 1)
    ctr += 1
ax_counts.set_ylim([-50, 800])
fig_counts.suptitle('Distribution of Number of Contigs within each cluster\n')
ax_counts.set_title('\n(The Red text indicates the total number of clusters)', fontsize = 20)
ax_counts_2.set_yticklabels([])
ax_counts_2.set_xlim([4.5, 28.5])
ax_counts.set_xlim([4.5, 28.5])

fig_counts.tight_layout()

fig_counts.savefig('Cluster_Analysis_Ward_Distribution_Contigs.pdf')


# In[ ]:


complete = df_cluster[df_cluster.columns[df_cluster.columns.str.startswith("Complete")]]
cols = complete.columns.tolist()
oplist = []
for c in cols:
    Cophenetic_Dist = float(c.replace("Complete_",""))
    df_grouped = complete[[c]].reset_index().groupby(c).count()
    singleton = len(df_grouped[df_grouped['ContigID'] <= 3])
    non_singleton = len(df_grouped[df_grouped['ContigID'] > 3])
    d = {'Cophenetic Distance':Cophenetic_Dist,
         '#Singleton':singleton,
         '#Non Singleton':non_singleton}
    oplist.append(d)
df_op = pd.DataFrame(data = oplist)
df_op = df_op.set_index('Cophenetic Distance')
df_op = df_op.sort_index()
indices = np.round(np.arange(0.4, 0.67, 0.01), 2)
df_op = df_op.loc[indices]


# In[ ]:


ax_singleton = df_op.plot(kind = 'bar', figsize = (16,10), stacked = True)
fig_singleton = ax_singleton.get_figure()
fig_singleton.tight_layout()
fig_singleton.savefig('Cluster_Singleton.pdf')


# In[ ]:





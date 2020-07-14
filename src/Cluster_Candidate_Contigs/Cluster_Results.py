import numpy as np
import pandas as pd
from os import listdir

wd = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/'
contigs = listdir(wd+'Mash/')
contigs.sort()
contigs = [c[:-6] for c in contigs]

cluster_files = listdir(wd+'Clustering_Results_Distance/')
cluster_id_files = []
for c in cluster_files:
        if ((c.startswith('Distance_')) and 
        	(c.endswith('.npy')) and 
        	('Matrix' not in c)):
                cluster_id_files.append(c)

df = pd.DataFrame(data = {'Contigs':contigs})

for c in cluster_id_files:
        c = c[:-4]
        split = c.split('_')
        d = round(float(split[2]),2)
        criterion = split[1]
        cluster_ids = np.load(wd+"Clustering_Results_Distance/"+c+".npy")
        df[criterion+'_'+str(d)] = cluster_ids
        print(criterion, d)

df = df.set_index('Contigs')
df.to_csv('/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Clustering_Results.csv')
print (df.max(axis = 0))



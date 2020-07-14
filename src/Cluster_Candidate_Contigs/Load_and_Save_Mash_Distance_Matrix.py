import pandas as pd
import numpy as np

filepath = '/fs/cbcb-scratch/hsmurali/Potential_Phages/Mash_k_11.txt'

df_dist = pd.read_csv(filepath, sep = ',', index_col = 'Reference', dtype = str)
del df_dist['Unnamed: 0']
ind = df_dist.columns.tolist()
ind.sort()
df_dist = df_dist.loc[ind, ind]

print (df_dist.info())

#Dist_Matrix = np.random.uniform(0,1,size = (78147, 78147))#df_dist.values

#np.save('/fs/cbcb-scratch/hsmurali/Potential_Phages/Random_Matrix',Dist_Matrix) 

import numpy as np
import pandas as pd
from os import listdir

df_distance = pd.DataFrame()

filepath = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash_Results_k_11_Dist/'

files = listdir(filepath)
files.sort()

Total_Number_Of_Files = 78938
op_list = np.zeros((Total_Number_Of_Files, Total_Number_Of_Files))
print("Matrix Created")
ctr = 0

columns = []

for f in files:
	print(ctr, f)
	df = pd.read_csv(filepath+f, sep = ',', index_col = 'Reference')
	del df['Unnamed: 0']
	indices = []

	if ctr == 0:
		columns = df.columns.tolist()
		columns.sort()
	df = df[columns]
	val = df.values

	assert len(columns) == Total_Number_Of_Files, "Assertion Failed "+str(len(columns)) 
	index = df.index.tolist()
	
	for ind in index:
		array_loc = columns.index(ind)
		indices.append(array_loc)
		#assert sum(op_list[array_loc]) == 0, "Assertion Failed.. Collision..." 
		#op_list[array_loc] = df.loc[ind].tolist()
	op_list[indices] = val
	del df
	ctr += 1

np.save('/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Distance_Matrix_K_11', op_list)



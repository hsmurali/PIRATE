import os
import sys
import numpy as np
import pandas as pd

Mash_Executable_Path = '/fs/cbcb-scratch/hsmurali/Mash-2.1.1/'
Data_Path = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash_Sketch/'
op_path = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash_Results_k_11_Dist/'

try:
	os.mkdir(op_path)
except Exception:
	pass

file1 = Data_Path + sys.argv[1]

files = os.listdir(Data_Path)

os.chdir(Mash_Executable_Path)

op_string = []

print (file1)

counter = 0

for f in files:
	counter += 1
	file2 = Data_Path + f
	command = './mash dist ' + file1 + ' '+file2
	results = os.popen(command).read()
	op_string += results.split('\n')

df = pd.DataFrame([x.split('\t') for x in op_string], columns = ['Reference','Query', 'Dist', 'P-Val', 'Matching-Kmer'])


df['Reference'] = df['Reference'].str.replace('/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash/', '')
df['Query'] = df['Query'].str.replace('/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash/', '')
df['Reference'] = df['Reference'].str.replace('.fasta','')
df['Query'] = df['Query'].str.replace('.fasta','')

df_pivoted = pd.pivot_table(df, values = 'Dist', index = ['Reference'], columns = ['Query'], aggfunc = np.sum)
df_pivoted = df_pivoted.reset_index()
df_pivoted.to_csv(op_path+sys.argv[1][:-4]+'.txt')

'''f_out = open(op_path + sys.argv[1][:-4] + '.txt', 'w')
f_out.write(op_string)
f_out.close()

compress_command = 'gzip ' + op_path + sys.argv[1][:-4]+'.txt'
results = os.popen(compress_command).read()'''


import numpy as np
import pandas as pd
from Bio import SeqIO
from os import mkdir


df = pd.read_csv('Filtered_Phage_List.csv')
print(len(df))
samples = np.unique(df['SampleID'].tolist())
print(len(samples))

sample_path = '/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/'
mash_path = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash/'

try:
	os.mkdir(mash_path)
except Exception:
	pass

for s in samples:
	break
	print(s)
	sample_id = s
	fastapath = sample_path + sample_id + '/' + sample_id + '.fna'
	record_dict = SeqIO.to_dict(SeqIO.parse(fastapath, "fasta"))	
	temp_df = df[df['SampleID'] == sample_id]
	contigs = temp_df['ContigID'].tolist()
	
	for contig in contigs:
		op_path = mash_path + sample_id + '_' + contig + '.fasta'
		SeqIO.write(record_dict[contig], op_path, "fasta")



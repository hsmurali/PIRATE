import sys
import subprocess
import numpy as np
import pandas as pd
import networkx as nx

from Bio import SeqIO
from packaging import version
from functools import partial
from os import listdir, mkdir
from os.path import isdir, isfile

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def cmd_exists(cmd):
	'''
	Function to check if a linux cmd exists
	Input:
		cmd: Linux Command
	Output:
		True if command exists, false otherwise. 
	'''
	return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def Check_Dependencies():
	'''
	Checks if all the dependencies are met! 
	'''
	pd_version = pd.__version__
	np_version = np.__version__
	nx_version = nx.__version__

	if version.parse(pd_version) < version.parse("1.2.4"):
		print("Incompatible Pandas version. Install Pandas (>=1.2.4)")
		return False
	if version.parse(np_version) < version.parse("1.18.3"):
		print("Incompatible Numpy version. Install Numpy (>=1.18.3)")
		return False
	if version.parse(nx_version) < version.parse("2.2"):
		print("Incompatible Networkx version. Install Networkx (>=2.2)")
		return False
	if not cmd_exists('blastn'):
		print("BLAST doesn't exist. Please install BLAST...")
		sys.exit(1)
	return True

def Run_BLAST(query_seq_path, num_threads, out_path):
	blast_cmd = ("blastn -query " + query_seq_path + " -db " +  "nt" + " -num_threads " + str(num_threads) + " -evalue 1e-3 "
				" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"" +
				" -out "+out_path)
	result = subprocess.Popen(blast_cmd,shell=True).wait()

def Write_Fasta_Seqs(d, out_path):
	out_buf = open(out_path,'w')
	for k in d:
		out_buf.write(">"+k+"\n")
		out_buf.write(d[k]+"\n")
	out_buf.close()
	
def Extract_Seqs(contig_path, repeats_path, prefix):
	d = {}
	fasta_sequences = SeqIO.parse(open(contig_path),'fasta')
	for f in fasta_sequences:
		d[f.name] = str(f.seq)

	repeats = open(repeats_path, 'r').readlines()
	out = {}
	for r in repeats:
		r = r.replace("\n","")
		out[prefix+r] = d[r]
	out

def Get_BLAST_Hits_Stats(grp):
	num_hits = len(grp)
	qlen = grp.iloc[0]['qlen']
	qstarts = grp['qstart'].tolist()
	qends = grp['qend'].tolist()
	cov = np.zeros(qlen)

	for i in range(0, len(qstarts)):
		start = min(qstarts[i], qends[i])
		end = max(qstarts[i], qends[i])
		cov[start:end] = 1
	return pd.Series({'#Hits':num_hits, 'Len':qlen, 'Breadth':cov.sum()*100.0/qlen})

def Summarize_BLAST_Results(blast_out_path):
	df = pd.read_csv(blast_out_path, sep = "\t", 
					 names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
							  'gapopen', 'qstart', 'qend', 'sstart', 'send',
							  'evalue', 'bitscore', 'qlen', 'slen'])
	df_grp = df.groupby('qseqid').apply(Get_BLAST_Hits_Stats)
	return df_grp

def Qualify_Alignments(row, contain, overhang_cutoff):
	if ((row['Overhang'] > overhang_cutoff) and (row['Overhang']/row['AlignLength'] > 0.1) and 
		(row['AlignLength']/min(row['QLen'], row['SLen']) < contain)):
		return "Internal_Match"
	
	elif ((row['QStart'] <= row['SStart']) and  
		  (row['AlignLength']/min(row['QLen'], row['SLen']) >= contain)):
		return "q in s"
	
	elif ((row['QStart'] >= row['SStart']) and 
		  (row['AlignLength']/min(row['QLen'], row['SLen']) >= contain)):
		return "s in q"
	
	elif row['QStart'] > row['SStart']:
		return "q overlaps s"
	
	else:
		return "s overlaps q"
	
def Build_Graph(df_Containments):
	edge_list = list(zip(df_Containments['Query'].tolist(), df_Containments['Subject'].tolist()))
	G = nx.DiGraph(edge_list)
	nodes = df_Containments['Query'].tolist()+df_Containments['Subject'].tolist()
   
	d_length = dict(zip(nodes, df_Containments['QLen'].tolist()+df_Containments['SLen'].tolist()))
	nx.set_node_attributes(G, d_length, name="length")

	d_edgetype = dict(zip(edge_list, df_Containments['Alignment_Type'].tolist()))
	nx.set_edge_attributes(G, d_edgetype, name="type")

	d_qcov = dict(zip(edge_list, df_Containments['QCov'].tolist()))
	nx.set_edge_attributes(G, d_qcov, name="qcov")

	d_scov = dict(zip(edge_list, df_Containments['SCov'].tolist()))
	nx.set_edge_attributes(G, d_scov, name="scov")

	d_qstart = dict(zip(edge_list, df_Containments['QStart'].tolist()))
	nx.set_edge_attributes(G, d_qstart, name="qstart")

	d_qend = dict(zip(edge_list, df_Containments['QEnd'].tolist()))
	nx.set_edge_attributes(G, d_qend, name="qend")

	d_sstart = dict(zip(edge_list, df_Containments['SStart'].tolist()))
	nx.set_edge_attributes(G, d_sstart, name="sstart")

	d_send = dict(zip(edge_list, df_Containments['SEnd'].tolist()))
	nx.set_edge_attributes(G, d_send, name="send")
	
	return G

def Simplify_Containment_Clusters(g, cluster_id):
	containment_clusters = []
	d_lengths = nx.get_node_attributes(g, "length")
	nodes, lengths = np.array(list(d_lengths.keys())), np.array(list(d_lengths.values()))
	nodes_sorted = nodes[np.argsort(lengths)[::-1]]
	
	visited = set({})
	
	for n in nodes_sorted:
		if n not in visited:
			cluster = [n] + list(g.neighbors(n))
			op = []
			for c in cluster[1:]:
				d = {'RepresentativeContig':n, 'Contig':c, 'EdgeType':g.edges[(n,c)]['type'], 
					 'RepresentativeLength':g.nodes[n]['length'], 'ContigLength':g.nodes[c]['length'],
					 'qstart':g.edges[(n,c)]['qstart'], 'qend':g.edges[(n,c)]['qend'], 
					 'sstart':g.edges[(n,c)]['sstart'], 'send':g.edges[(n,c)]['send'],
					 'Cluster_ID':cluster_id}
				op.append(d)
			cluster_id += 1
			g.remove_nodes_from(cluster)
			visited.update(cluster)
			if len(cluster) > 1:
				containment_clusters+= op
	return containment_clusters

def Load_PAF(filepath, percent_identity_cutoff, containment_cutoff, overhang_cutoff):
	header = ['Query','QLen','QStart','QEnd','Orientation','Subject',
			  'SLen','SStart','SEnd','Matches','AlignLength','MAPQ']
	op = []
	with open(filepath) as fileobject:
		for l in fileobject:
			l = l.split('\t')[:12]
			op.append(dict(zip(header, l)))
	df = pd.DataFrame(op)
	df[['QLen','QStart','QEnd','SLen','SStart',
		'SEnd','Matches','AlignLength','MAPQ']] = df[['QLen','QStart','QEnd','SLen','SStart',
													  'SEnd','Matches','AlignLength','MAPQ']].astype('int')
	df['PIdent'] = df['Matches']/df['AlignLength']*100
	df['QCov'] = df['AlignLength']/df['QLen']*100
	df['SCov'] = df['AlignLength']/df['SLen']*100
	df['Overhang'] = (np.minimum(df['SStart'],df['QStart']) + np.minimum(df['QLen']-df['QEnd'], 
																		 df['SLen']-df['SEnd']))
	df = df[df['Query'] != df['Subject']]
	df_filtered = df.loc[df['PIdent'] >= percent_identity_cutoff, :].copy()
	df_filtered = df_filtered.loc[df_filtered.groupby(['Query','Subject'])['AlignLength'].idxmax()]
	df_filtered.loc[:,'Alignment_Type'] = df_filtered.apply(partial(Qualify_Alignments, 
																	contain = containment_cutoff,
																	overhang_cutoff = overhang_cutoff), axis = 1)

	return df_filtered

def Filter_Seq_Description(desc, keywords, delim):
	desc = desc.lower()
	desc_arr = desc.split(delim)
	for d in desc_arr:
		for k in keywords:
			if k in d:
				return 1
	return 0

def Get_Contig(name):
	x = name.split('_')
	return "_".join(x[:-1])

def Get_Sample(name):
	return name.split('.')[0]

def Plot_Genes(Annotations, title):
	colors = ['red', 'blue', 'green', 'grey', 'orange', 'gold', 'purple', 'olive', 
			  'cyan', 'magenta', 'lime', 'teal', 'black', 'yellow', 'lightcoral', 'brown', 
			  'peachpuff', 'plum', 'peru', 'sienna', 'maroon', 'indigo', 'darkolivegreen', 
			  'crimson', 'thistle', 'dodgerblue', 'firebrick', 'mediumslateblue', 'lightsalmon', 
			  'khaki', 'darkseagreen', 'lawngreen', 'sandybrown', 'tan', 'royalblue', 'goldenrod', 
			  'springgreen', 'cadetblue', 'pink', 'slategrey']
	Annotations = Annotations.sort_values(by = ['sstart','ContigLength'], ascending = [True,False])
	contigs = Annotations['Contig'].tolist()+Annotations['RepresentativeContig'].tolist()
	contig_starts = Annotations['sstart'].tolist()
	contig_ends = Annotations['send'].tolist()
	rep_len = int(Annotations.iloc[0]['RepresentativeLength'])

	preferred_names = Annotations['Preferred_name'].tolist()
	full_annotations = Annotations['Description'].tolist()

	unique_annotations = []
	for i in range(0, len(full_annotations)):
		if len(full_annotations[i].split(" ")) >= 20:
			unique_annotations.append(preferred_names[i])
		else:
			unique_annotations.append(full_annotations[i])
	unique_annotations = np.unique(unique_annotations)
	Annotation_keys = dict(zip(unique_annotations,
							   [np.zeros(rep_len) for i in range(len(unique_annotations))]))
	
	fig = plt.figure(figsize=(16, 8))
	gs = GridSpec(nrows=2, ncols=1, height_ratios=[3.5, 1])
	ax0 = fig.add_subplot(gs[1, 0])
	ax1 = fig.add_subplot(gs[0, 0])
	
	for i in range(len(Annotations)):
		try:
			start = int(Annotations.iloc[i]['e_qstart']+Annotations.iloc[i]['sstart'])
			end = int(Annotations.iloc[i]['e_qend']+Annotations.iloc[i]['sstart'])
			annotation = Annotations.iloc[i]['Description']
			if len(annotation.split(" ")) >= 20:
				annotation = Annotations.iloc[i]['Preferred_name']
			Annotation_keys[annotation][start:end] = 1
		except ValueError:
			continue
			
	for x in range(len(Annotation_keys)):
		if str(unique_annotations[x]) == '-':
			continue
		X = Annotation_keys[unique_annotations[x]]
		X[X == 0] = np.inf
		X[X == 1] = x
		ax0.plot(X, color = colors[(x+5)%len(colors)], linewidth = 6, label = unique_annotations[x])
		
	ax0.plot([0, rep_len],[len(unique_annotations)]*2, color = 'black', linewidth = 10)
	ax0.set_yticks([])
	
	d_samples = {}
	max_so_far = 0
	for i in range(len(contig_starts)):
		sample = contigs[i].split('_')[0]
		try:
			ctr = d_samples[sample]
		except:
			d_samples[sample] = max_so_far
			ctr = max_so_far
			max_so_far += 1
		ax1.plot([contig_starts[i], contig_ends[i]], [ctr, ctr], color = 'grey')
	ax1.plot([0, rep_len],[len(d_samples)+2]*2, color = 'black', linewidth = 10)
	ax1.set_yticks([])
	ax1.set_xticks([])
	fig.tight_layout()
	ax1.set_title(title)
	
	handles, labels = ax0.get_legend_handles_labels()
	lgd = fig.legend(handles, labels, loc=1, ncol = 2, fontsize = 17, 
			   bbox_to_anchor = (-0.125, -1.0, 1, 1), frameon = False)
	
	return fig, lgd
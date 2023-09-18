import argparse as ap
from PIRATE_Utils import *

if __name__ == "__main__":
	parser = ap.ArgumentParser(description="PIRATE: Phage Identification fRom Assembly-graph varianT Elements."+
											"PIRATE is a pipeline for identifying phage elements from assembly graph variants such"+ 
											" as interspersed repeats and indel bubblesWe assume MetaCarvel was run on metagenomic assemblies "+
											"and variants were called out. After filtering out contigs found in databases with sufficient coverage, we run PIRATE.")

	parser.add_argument("-i", "--input", help="A fasta file of sequences to cluster", required=True)
	parser.add_argument("-d", "--percent_identity_cutoff", help="percent identity to keep contigs within the same cluster.(default = 95)", default = 95, required=False)
	parser.add_argument("-c", "--containment_cutoff", 
						help="Fraction of the query contig to be covered by the centroid contig to retain within the cluster.(default = 0.80)", 
						default = 0.80, required=False)
	parser.add_argument("-o", "--overhang", 
						help="Length of the query contig overhanging with respect to the cluster centorid. (default = 100bp)", 
						default = 100, required=False)
	parser.add_argument("-t","--num_threads", help="Number of threads. (default = 8)", default = 8, required=False)
	parser.add_argument("-out","--output_directory", help="Location to write the outputs to", required = True)
	parser.add_argument("-pre","--prefix", help = "Prefix for output files. (default = PIRATE)", default = "PIRATE", required = False)
	args = parser.parse_args()

	Check_Dependencies()
	
	if not isfile(args.input):
		print("File ",args.input," not found... Check file path...")

	if not isdir(args.output_directory):
		mkdir(args.output_directory)

	percent_identity_cutoff = float(args.percent_identity_cutoff)
	containment_cutoff = float(args.containment_cutoff)
	overhang_cutoff = float(args.overhang)
	threads = int(args.num_threads)
	command_all_vs_all = "minimap2 -X -w 5 -t "+str(threads)+" "+ args.input+" "+args.input+" > "+args.output_directory+"/"+args.prefix+".all-vs-all.paf"
	result = 	subprocess.Popen(command_all_vs_all,shell=True).wait()

	df = Load_PAF(args.output_directory+"/"+args.prefix+".all-vs-all.paf", percent_identity_cutoff, containment_cutoff, overhang_cutoff)
	df_Containments = df.loc[(df['Alignment_Type'] == 'q in s')|(df['Alignment_Type'] == 's in q'), :].copy()
	G = Build_Graph(df_Containments)
	
	conn = list(nx.weakly_connected_components(G))
	containment_clusters = []
	cluster_id = 0

	for c in conn:
		g = nx.Graph(G.subgraph(c))
		clusters = Simplify_Containment_Clusters(g, cluster_id)
		if len(clusters) > 1: cluster_id = clusters[-1]['Cluster_ID']+1
		containment_clusters += clusters
	df_containment_clusters = pd.DataFrame(containment_clusters)
	df_containment_clusters.to_csv(args.output_directory+"/"+args.prefix+".clusters", sep = "\t")

	nx.write_gml(G, args.output_directory+"/"+args.prefix+".containment.gml")
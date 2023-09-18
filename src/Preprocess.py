import argparse as ap
from PIRATE_Utils import *

if __name__ == "__main__":
	parser = ap.ArgumentParser(description="PIRATE: Phage Identification fRom Assembly-graph varianT Elements."+
											"This is the pieline to filter potential repeat contigs of interest that is "+
											"likely to contain novel phage-elements. To that end, we align all contigs identified as "+
											"interspersed repeats against the non-redundant -nt database and retain only "+
											"poorly covered contigs for further downstream tasks. ")
	parser.add_argument("-i", "--input", help="A fasta file of assembled contigs", required=True)
	parser.add_argument("-l", "--contig_list", help="A test file of contigs to be considered for putative phage discovery", required=True)
	parser.add_argument("-o", "--output_directory", help="Directory to write outputs to.", required=True)

	parser.add_argument("-c", "--breadth_coverage_cutoff", help="breadth of coverage to filter contigs. (default = 50)", default = 50, required=False)
	parser.add_argument("-len", "--min_len", help="minimum length of contigs to considered to retain for further analyses. (default = 3000)", default = 3000, required=False)
	parser.add_argument("-t","--num_threads", help="Number of threads. (default = 8)", default = 8, required=False)
	
	parser.add_argument("-pre","--prefix", help = "Prefix for output files (default = PIRATE.preprocess)", default = "PIRATE.preprocess", required = False)
	
	args = parser.parse_args()

	Check_Dependencies()
	
	if not isfile(args.input):
		print("File ",args.input," not found... Check file path...")
		sys.exit(1)

	if not isfile(args.contig_list):
		print("File ",args.contig_list," not found... Check file path...")
		sys.exit(1)
	
	if not isdir(args.output_directory):
		mkdir(args.output_directory)

	breadth_coverage_cutoff = float(args.breadth_coverage_cutoff)
	min_len = float(args.min_len)
	threads = int(args.num_threads)
	prefix = args.output_directory+args.prefix
	
	out = Extract_Seqs(args.input, args.contig_list, args.prefix+".")
	Write_Fasta_Seqs(out, prefix+'.repeats.fa')
	Run_BLAST(prefix+'.repeats.fa', threads, prefix+'.BLAST.out')
	cmd = 'gzip '+prefix+'.BLAST.out'
	result = subprocess.Popen(cmd,shell=True).wait()

	df_blast_summary = Summarize_BLAST_Results(prefix+'.BLAST.out.gz')
	df_blast_summary_filtered = df_blast_summary[(df_blast_summary['Len'] > min_len) & 
												 (df_blast_summary['Breadth'] <= breadth_coverage_cutoff)]
	
	print(args.prefix, len(df_blast_summary_filtered))
	df_blast_summary.to_csv(prefix+'.BLAST.summary', sep = "\t")
	df_blast_summary_filtered.to_csv(prefix+'.BLAST.filtered.summary', sep = "\t")
	
	putative_contigs = df_blast_summary_filtered.index.tolist()
	o = open(prefix+'.PIRATE.contigs.list','w')
	for p in putative_contigs:
		o.write(p+"\n")
	o.close()

	d_repeats_pirate = Extract_Seqs(prefix+'.repeats.fa', prefix+'.PIRATE.contigs.list', "")
	Write_Fasta_Seqs(d_repeats_pirate, prefix+'.PIRATE.contigs.fa')

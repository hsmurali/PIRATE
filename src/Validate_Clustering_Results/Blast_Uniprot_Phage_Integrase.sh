#!/bin/bash
#SBATCH -J BLAST_UniProt # Job name
#SBATCH -o BLAST_UniProt.o # Name of output file
#SBATCH -e BLAST_UniProt.e # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb

module load blast
source activate hsmurali27
makeblastdb -in /fs/cbcb-scratch/hsmurali/Phage_Discovery_Project/uniprot-organism_phage+integrase.fasta  -parse_seqids -dbtype prot
blastx -query /fs/cbcb-scratch/hsmurali/Phage_Discovery_Project/Phage_Candidate_Contigs_New.fasta -db /fs/cbcb-scratch/hsmurali/Phage_Discovery_Project/uniprot-organism_phage+integrase.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -out /fs/cbcb-scratch/hsmurali/Phage_Discovery_Project/Phage_Candidates_Genes_New.txt
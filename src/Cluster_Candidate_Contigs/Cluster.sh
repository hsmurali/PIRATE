#!/bin/bash
#SBATCH -J Cluster_Contigs # Job name
#SBATCH -o Cluster_Contigs.o # Name of output file
#SBATCH -e Cluster_Contigs.e # Name of error file
#SBATCH --mail-user=hsmurali@cs.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=5-18:00:00
#SBATCH --qos=xlarge
#SBATCH --mem=450gb

source activate hsmurali27
python2 Perform_Clustering.py




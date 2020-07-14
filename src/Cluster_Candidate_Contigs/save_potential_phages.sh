#!/bin/bash
#SBATCH -J id_contigs # Job name
#SBATCH -o id_contigs.o%j # Name of output file
#SBATCH -e id_contigs.e%j # Name of error file
#SBATCH --mail-user=hsmurali@cs.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=5
#SBATCH --qos=large
#SBATCH --mem=128gb

source activate hsmurali27
python2 /fs/cbcb-scratch/hsmurali/Potential_Phages/save_potential_phages_fasta.py


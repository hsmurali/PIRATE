#!/bin/bash
#SBATCH -J Mash # Job name
#SBATCH -o Mash.o%j # Name of output file
#SBATCH -e Mash.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=1-790

source activate hsmurali27
file=`head -n ${SLURM_ARRAY_TASK_ID} Mash_Commands.txt | tail -n 1`
python2 $file

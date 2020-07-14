#!/bin/bash
#SBATCH -J Build_Dist_Matrix # Job name
#SBATCH -o Build_Dist_Mat.o # Name of output file
#SBATCH -e Build_Dist_Mat.e # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=2-18:00:00
#SBATCH --qos=xlarge
#SBATCH --mem=450gb

source activate hsmurali27
python2 Build_Distance_Matrix.py

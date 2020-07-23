
<h1>Clustering Candidate Phage Contigs using MASH</h1>

In this phase of the pipeline we cluster the contigs identified in the previous step to identify phage elements. This part of the code assumes the availability of a cluster. run the following programs in the following order. This pipeline is still a work in progress and we hope to clean this up soon, and tie together a single script that gets the job done. the python scripts are mostly a wrapper over the MASH software. 

1. **save_potential_phages.sh** - This script creates a directory "Mash" and creates a fasta file for each of the candidate contigs. This script invokes the python script "save_potential_phages_fasta.py".


2. **make_mash_sketches.py** - This script sketches from the fasta files generated in the previous step. The script runs pretty fast and creates a directory "Mash_Sketch" and dumps the sketch files here. 


3. **Prepare_Mash_Sketch_Distance_Commands.py** - This python script generates "Mash_Commands.txt" which repeatedly invokes the "run_mash_sketch_distance.py". 


4. **Submit_Mash.sh** - This bash script invokes "run_mash_sketch_distance.py" using "Mash_Commands.txt" and places the pairwise mash distance calculated in directory "Mash_Results_K_11_Dist".


5. **Build_Distance_Matrix.sh** - This bash script invokes "Build_Distance_Matrix.py" which creates the pairwise-mash distance matrix, which is then used to clusters the contigs. 


6. **Cluster.sh** - This script invokes the script Perform_Clustering.py which clusters the contigs based on different clustering methods namely single linkage, complete linkage, average linkage, Ward's, centroid and weighted linkage. 


7. **Cluster_Results.py** - This script smmarizes the output from the previous step and writes output to a ".csv" file which is used for validation. 

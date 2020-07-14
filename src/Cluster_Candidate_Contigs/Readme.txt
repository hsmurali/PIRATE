Please run the following commands in order:


1. Build Mash Sketches(each mash sketch contains sketch of 100 fasta files, assuming each fasta file contains one contig)
		python make_mash_sketches.py (you may just submit it as any job in SLURM, this runs very fast)
			Writes output into Mash_Sketch folder. 

2. Create a text file containing a list of commands to run compute mash distance. This step is very important to submit an array job on the SLURM Supercomputer.
   Invokes the python script run_mash_sketch_distance script repeatedly from SLURM to compute the distance values. 
		python Prepare_Mash_Sketch_Distance_Commands.py
			Writes output to Mash_Commands.txt

3. Perform distance computation in the clsuter as an array job. Loads commands from the output of last step, i.e. Mash_Commands.txt
		sbatch Submit_Mash.sh
			Writes output as compressed text files into Mash_Results_k_11 directory

4. Prepare commands to process the compressed text output files from previous step as an array job. This step is very important to submit an array job on the SLURM Supercomputer.
   Invokes Merge_Mash_Op_Large_Batch.py repeatedly.
		python Prepare_Dist_Commands.py
			Writes output into Mash_Commands)Extract_Dist.txt

5. Submit a slurm job to extract distance information from the compressed text files from step 3. 
		sbatch Merge_Mash_Large.sh
			Writes output into Mash_Results_k_11_Dist directory.

6. Build the symmetric Distance Matrix that can be used to perform clustering. Invokes Build_Distance_Matrix.py.
		sbatch Build_Distance_Matrix.sh
			Writes output into a ".npy" file, Distance_Matrix_k_11.npy

7. Perform Clustering with the above distance matrix. Invokes the Perform_Clustering.py script.
		sbatch Cluster.sh
			Writes output into the working directory for different linkage methods and cophenetic-distance cutoff values. 

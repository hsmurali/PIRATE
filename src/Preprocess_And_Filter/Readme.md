
<h1> Preprocess and Filter Data </h1>

This step of the pipeline PIRATE is aimed at preprocessing the blast hits against the non-redundant (NR) database on all the contigs. The script assumes that metacarvel is already run on all the samples in which you would want to identify Phage/prophage elements. To learn how to run metacarvel please check out, https://github.com/marbl/MetaCarvel. After running metacarvel BLAST the contigs against the NR database and then run the following scripts in order for an exploratory analysis of the same. 

1. **Load_and_Pickle_BLAST_Outputs.py**- This script is aimed at loading the BLAST outputs and save them as Pickle objects for easier loading. 



2. **Compute_Breadth_of_Coverage.py**- This script computes the breadth of coverage. We use this to primariy identify contigs that are outliers. 



3. **Analyze_BLAST_Hits.py** - This is an optional script aimed at preparing plots that shows the various BLAST parameters. 



4. **Preprocess_KRAKEN_Outputs.py** -  We also ran KRAKEN on the contigs to compute taxa labels on the contigs and removed contigs that had a taxa label below family. This did not influence the accuracy of the phage candidate contigs. 



5. **Identify_Candidate_Phage_Contigs.py** - Based on the coverage scores we have computed in the previous step this script identifies a list of contigs considered for clustering in the subsequent steps. 


```python

```

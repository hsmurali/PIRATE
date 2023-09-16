#!/bin/bash

input_seqs=${1}
output_dir=${2}
prefix=${3}

mkdir ${output_dir}
mkdir ${output_dir}/Prodigal/
prodigal -i ${input_seqs} \
         -p meta \
         -a ${output_dir}/Prodigal/${prefix}.prodigal.faa \
         -d ${output_dir}/Prodigal/${prefix}.prodigal.fna \
         -o ${output_dir}/Prodigal/${prefix}.prodigal.out

eggnog_db_path=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/eggnog_functional_comparison/databases/
mkdir ${output_dir}/EggNOG/

python3 -c "import sqlite3; print(sqlite3. sqlite_version)"
emapper.py --version --data_dir ${eggnog_db_path}
emapper.py -m diamond --no_annot --no_file_comments --cpu 24 \
		   --data_dir ${eggnog_db_path} \
           -i ${output_dir}/Prodigal/${prefix}.prodigal.faa \
           --output ${output_dir}/EggNOG/${prefix}.eggnog.out

emapper.py -m no_search \
		   --annotate_hits_table ${output_dir}/EggNOG/${prefix}.eggnog.out.emapper.seed_orthologs \
           --no_file_comments -o ${output_dir}/EggNOG/${prefix}.eggnog.out \
           --cpu 24 \
           --data_dir ${eggnog_db_path} --dbmem
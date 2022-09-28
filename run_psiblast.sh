#!/bin/sh
PGPATH=/Users/suhancho/program/external/ncbi-blast-2.13.0+/bin/
FILEPATH=/Users/suhancho/data/Uniprot_metalbinding_challenge/
# DB_POS=${FILEPATH}/ref/uniref50.fasta
DB_POS=/Users/suhancho/data/Uniprot_metalbinding_challenge//ref/random_split_uniref/uniref.10percent.1.fasta

mkdir -p ${FILEPATH}/PSSM/
while read fasta_file
do
protname=$(echo ${fasta_file} | cut -f7 -d/)
${PGPATH}/psiblast \
-query ${fasta_file} \
-db ${DB_POS} \
-num_iterations 2 \
-evalue 0.001 \
-num_threads 48 \
-out_ascii_pssm ${FILEPATH}/PSSM/${protname}.txt
done < ${FILEPATH}/pos_list.txt
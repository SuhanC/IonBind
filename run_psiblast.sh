#!/bin/sh
PGPATH=/Users/suhancho/program/external/ncbi-blast-2.13.0+/bin/
FILEPATH=/Users/suhancho/data/Uniprot_metalbinding_challenge/
DB_POS=${FILEPATH}/POS_TRAIN_FULL.fasta.pto
mkdir -p ${FILEPATH}/PSSM/

while read fasta_file
do
${PGPATH}/psiblast \
-query fasta_file \
-db ${DB_POS} \
-num_iterations 2 \
-evalue 0.001 \
-num_threads 8 \
-out_pssm ${FILEPATH}/PSSM/fasta_file.txt

done < ${FILEPATH}/pos_list.txt
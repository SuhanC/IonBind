#!/bin/sh
INFILE=/Users/suhancho/data/Uniprot_metalbinding_challenge/ref/uniref50.fasta
OUTFILE=/Users/suhancho/data/Uniprot_metalbinding_challenge/ref/random_split_uniref/uniref.10percent

for i in 1 2 3 4 5 6 7 8 9 10 
do
pyfastx sample \
-o ${OUTFILE}.${i}.fasta -p 0.1 ${INFILE}
done

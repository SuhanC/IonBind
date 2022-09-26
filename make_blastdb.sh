PGPATH=/Users/suhancho/program/external/ncbi-blast-2.13.0+/bin/
FASTA=/Users/suhancho/data/Uniprot_metalbinding_challenge/POS_TRAIN_FULL.fasta
${PGPATH}/makeblastdb -in ${FASTA} -parse_seqids -title "Cookbook demo" -dbtype prot

PGPATH=/Users/suhancho/program/external/ncbi-blast-2.13.0+/bin/
# FASTA=/Users/suhancho/data/Uniprot_metalbinding_challenge/ref/uniref50.fasta
FASTA=/Users/suhancho/data/Uniprot_metalbinding_challenge//ref/random_split_uniref/uniref.10percent.1.fasta

${PGPATH}/makeblastdb -in ${FASTA} -parse_seqids -title "uniref50_downloaded.split1" -dbtype prot

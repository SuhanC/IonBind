#!/bin/sh
PGPATH=/Users/suhancho/program/external/mmseqs/bin/
FILEPATH=/Users/suhancho/data/Uniprot_metalbinding_challenge/
DB_POS=${FILEPATH}/ref/uniref50.fasta
mkdir -p ${FILEPATH}/PSSM/
mkdir -p ${FILEPATH}/mmseq2_tmp/
mkdir -p ${FILEPATH}/mmseq2_out/
mkdir -p ${FILEPATH}/mmseq2_profile_out/
mkdir -p ${FILEPATH}/mmseq2_PSSM/
mkdir -p ${FILEPATH}/mmseq2_pos_sequence/

# ${PGPATH}/mmseqs createdb ${DB_POS} ${FILEPATH}/ref/uniref50.mmseq2.fasta 

DB_MMSEQS=${FILEPATH}/ref/uniref50.mmseq2.fasta
PSSMOUT=${FILEPATH}/mmseq2_PSSM/
MMSEQOUT=${FILEPATH}/mmseq2_out/
PROFILEOUT=${FILEPATH}/mmseq2_profile_out/
MMSEQTMP=${FILEPATH}/mmseq2_tmp/
MMSEQINPUT=${FILEPATH}/mmseq2_pos_sequence/

while read fasta_file
do
protname=$(echo ${fasta_file} | cut -f7 -d/)
${PGPATH}/mmseqs createdb ${fasta_file} ${MMSEQINPUT}/${protname}
${PGPATH}/mmseqs search ${MMSEQINPUT}/${protname} ${DB_MMSEQS} ${MMSEQOUT}/${protname}.out ${MMSEQTMP} -e 0.0005 -a --num-iterations 2
${PGPATH}/mmseqs result2profile ${MMSEQINPUT}/${protname} ${DB_MMSEQS} ${MMSEQOUT}/${protname}.out ${PROFILEOUT}/${protname}.profile
${PGPATH}/mmseqs profile2pssm ${PROFILEOUT}/${protname}.profile ${PSSMOUT}/${protname}.pssm
done < ${FILEPATH}/pos_list.txt
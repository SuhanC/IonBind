OUTPATH=/Users/suhancho/data/Uniprot_metalbinding_challenge/pos_sequence_single_fasta/
mkdir -p ${OUTPATH}

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outname=${line#>}
        outname_fasta=$(echo ${outname} | cut -f2 -d"|")
        outfile=${OUTPATH}/${outname_fasta}.fasta

        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < /Users/suhancho/data/Uniprot_metalbinding_challenge/POS_TRAIN_FULL.fasta
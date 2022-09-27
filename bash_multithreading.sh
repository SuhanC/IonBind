# #!/bin/sh
# th="16"                                             # 동시 생성할 최대 쓰레드 숫자 지정 #
# list=( 10 5 8 2 3 9 1 2 5 3 9 2 5 ) 
# for a in ${list[@]}
#   do j_count=`jobs -l|wc|awk '{print $1}'`
#   if [[ $j_count -ge $th ]];then
#     until [[ $j_count -lt $th ]]
#       do j_count=`jobs -l|wc|awk '{print $1}'`
#       sleep 0.1
#       done
#   fi
# ###########################################
# # threads 에 명령 넘기는 부분.
#   sleep $a &
# ###########################################
#   done
 
# lastPIDs=`jobs -l|awk '{print $2}'`
# wait $lastPIDs
# echo "";echo "work complet."


PGPATH=/Users/suhancho/program/external/ncbi-blast-2.13.0+/bin/
FILEPATH=/Users/suhancho/data/Uniprot_metalbinding_challenge/
DB_POS=${FILEPATH}/ref/uniref50.fasta
mkdir -p ${FILEPATH}/PSSM/




th="16"                                             # 동시 생성할 최대 쓰레드 숫자 지정 #
while read fasta_file
  do j_count=`jobs -l|wc|awk '{print $1}'`
  if [[ $j_count -ge $th ]];then
    until [[ $j_count -lt $th ]]
      do j_count=`jobs -l|wc|awk '{print $1}'`
      sleep 0.1
      done 
  fi
###########################################
# threads 에 명령 넘기는 부분.
  sleep 
  protname=$(echo ${fasta_file} | cut -f7 -d/)
  ${PGPATH}/psiblast \
  -query ${fasta_file} \
  -db ${DB_POS} \
  -num_iterations 2 \
  -evalue 0.001 \
  -num_threads 8 \
  -out_ascii_pssm ${FILEPATH}/PSSM/${protname}.txt &
###########################################
  done <${FILEPATH}/pos_list.txt
 
lastPIDs=`jobs -l|awk '{print $2}'`
wait $lastPIDs
echo "";echo "work complet."

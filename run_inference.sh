#!/bin/bash
END=1194
scriptpath=/home/suhan/script/IonBind/inference_runcode/

source /home/suhan/miniconda3/etc/profile.d/conda.sh
conda activate keras

for i in $(seq 0 $END);
do

python ${scriptpath}/job.${i}.py

echo job ${i} running

done
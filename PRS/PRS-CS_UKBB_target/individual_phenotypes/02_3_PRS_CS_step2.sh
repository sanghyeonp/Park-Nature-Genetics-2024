#!/bin/bash
source 00_global_variables.sh

OUT=PRS

for idx in "${!Phenotypes[@]}"
do
    outDIR="${OUT}/${Phenotypes[$idx]}"
    cd ${outDIR}
    echo ${outDIR}
    cat ./out/*.txt > ./score_sum.txt
done

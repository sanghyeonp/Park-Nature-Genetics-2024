#!/bin/bash
source 00_global_variables.sh

OUT=PRS
UKB_bim=ukb_eur_unrel_comb
for idx in "${!Phenotypes[@]}"
do
    outDIR="${OUT}/${Phenotypes[$idx]}"
    cd ${outDIR}
    plink --bfile ${UKB_bim} --score ./score_sum.txt 2 4 6 sum --out ./score_sum_indiv
done

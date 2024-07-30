#!/bin/bash
source 00_global_variables.sh

OUT=PRS
KoGES_bim=KCHIP.snp.imputed.HRC.rsID.rmdup
for idx in "${!Phenotypes[@]}"
do
    outDIR="${OUT}/${Phenotypes[$idx]}"
    cd ${outDIR}
    plink --bfile ${KoGES_bim} --score ./score_sum.txt 2 4 6 sum --out ./score_sum_indiv
done

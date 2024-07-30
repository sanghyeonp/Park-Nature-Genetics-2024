#!/bin/bash

trait=WC
n_thread=5
###

exe_gcta=gcta-1.94.1
bfile=UKB_random10k_noqc
file_leadsnp=../leadSNP/leadSNP_list.${trait}.txt
dir_cojo_input=../cojo_input

dir_out=cojo_out.${trait}
mkdir -p ${dir_out}

${exe_gcta} --bfile ${bfile} \
            --maf 0.01 \
            --diff-freq 0.2 \
            --cojo-file ${dir_cojo_input}/cojo_input.${trait}.txt \
            --cojo-slct \
            --cojo-cond ${file_leadsnp} \
            --out ${dir_out}/cojo.${trait} \
            --threads ${n_thread}
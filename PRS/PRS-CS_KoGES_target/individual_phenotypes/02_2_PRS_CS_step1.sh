#!/bin/bash
source 00_global_variables.sh

N_THREADS=16
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

OUT=./PRS
ref_1kg_eur=../ldblk_1kg_eur
KoGES_bim=KCHIP.snp.imputed.HRC.rsID.rmdup

for idx in "${!Phenotypes[@]}"
do
    OUT2="${OUT}/${Phenotypes[$idx]}"
    mkdir -p $OUT2

    OUT3="${OUT2}/out"
    mkdir -p $OUT3

    sumstat=${Sumstat_for_PRS_noUKB[$idx]}
    samplesize=${SampleSize_noUKB[$idx]}

    python3 ../PRScs/PRScs.py --ref_dir "${ref_1kg_eur}" \
                        --bim_prefix "${KoGES_bim}" \
                        --sst_file "${sumstat}" \
                        --n_gwas ${samplesize} \
                        --phi 1e-2 \
                        --out_dir "${OUT3}/${Phenotypes[$idx]}"
done

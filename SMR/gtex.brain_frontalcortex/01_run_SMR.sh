#!/bin/bash

gwas=MetS_GWAS.txt
ref_random_ukb=.test.rand10k
eqtl_besd=./Brain_Frontal_Cortex_BA9.lite

smr_Linux --bfile ${ref_random_ukb} \
    --gwas-summary ${gwas} \
    --beqtl-summary ${eqtl_besd} \
    --out ./smr.MetS.GTEx.Brain_Frontal_Cortex_BA9 \
    --smr-multi \
    --set-wind 1000 \
    --ld-multi-snp 0.2 \
    --thread-num 5

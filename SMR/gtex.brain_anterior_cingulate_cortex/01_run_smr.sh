#!/bin/bash

gwas=MetS_GWAS.txt
ref_random_ukb=.test.rand10k
eqtl_besd=./Brain_Anterior_cingulate_cortex_BA24/Brain_Anterior_cingulate_cortex_BA24

smr_Linux --bfile ${ref_random_ukb} \
    --gwas-summary ${gwas} \
    --beqtl-summary ${eqtl_besd} \
    --out ./smr.MetS.GTEx.Brain_Anterior_cingulate_cortex_BA24 \
    --smr-multi \
    --set-wind 1000 \
    --ld-multi-snp 0.2 \
    --thread-num 5


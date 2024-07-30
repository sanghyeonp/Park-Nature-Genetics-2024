#!/bin/bash

gwas=MetS_GWAS.txt
ref_random_ukb=.test.rand10k
eqtl_besd=./Muscle_Skeletal/Muscle_Skeletal

smr_Linux --bfile ${ref_random_ukb} \
    --gwas-summary ${gwas} \
    --beqtl-summary ${eqtl_besd} \
    --out ./smr.MetS.GTEx.Muscle_Skeletal \
    --smr-multi \
    --set-wind 1000 \
    --ld-multi-snp 0.2 \
    --thread-num 5


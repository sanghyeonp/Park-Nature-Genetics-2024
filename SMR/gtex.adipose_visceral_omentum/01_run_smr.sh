#!/bin/bash

gwas=MetS_GWAS.txt
ref_random_ukb=.test.rand10k
eqtl_besd=./Adipose_Visceral_Omentum/Adipose_Visceral_Omentum

smr_Linux --bfile ${ref_random_ukb} \
    --gwas-summary ${gwas} \
    --beqtl-summary ${eqtl_besd} \
    --out ./smr.MetS.GTEx.Adipose_Visceral_Omentum \
    --smr-multi \
    --set-wind 1000 \
    --ld-multi-snp 0.2 \
    --thread-num 5


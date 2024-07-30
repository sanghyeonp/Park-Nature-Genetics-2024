#!/bin/bash

gwas=MetS_GWAS.txt
ref_random_ukb=.test.rand10k
eqtl_besd=./BrainMetaV2/BrainMeta_cis_eQTL_chr

for chr in {1..22};do
    smr_Linux --bfile ${ref_random_ukb} \
        --gwas-summary ${gwas} \
        --beqtl-summary ${eqtl_besd}${chr} \
        --out ./smr.MetS.BrainMetaV2 \
        --smr-multi \
        --set-wind 1000 \
        --ld-multi-snp 0.2 \
        --thread-num 5
done
#!/bin/bash

###
trait=WC
target_cohort=UKB


###
mkdir -p ${trait}

dir_prsice=../PRSice
gwas=../gwas.${trait}_noUKB.PRSice_input.txt
target=./ukb_eur_unrel_comb
ld_ref=./EUR
pheno_file=./ukb.pheno.txt
cov_file=./ukb.cov.txt
mets_prev=0.253

Rscript ${dir_prsice}/PRSice.R \
    --out ./${trait}/${trait}.${target_cohort} \
    --prsice ${dir_prsice}/PRSice_linux \
    --base ${gwas} \
    --base-maf MAF:0.01 \
    --stat BETA \
    --beta \
    --ld ${ld_ref} \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --target ${target} \
    --pheno ${pheno_file} \
    --binary-target T \
    --prevalence ${mets_prev} \
    --cov ${cov_file} \
    --quantile 10 \
    --thread 10
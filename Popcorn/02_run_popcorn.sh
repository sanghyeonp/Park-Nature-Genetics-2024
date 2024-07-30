#!/bin/bash

eur=MetS_ukbincluded
eas=KoGES

s1=QC.MetS12_GWAS_MetS_maf0.01_popcornIN.tsv
s2=KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid.modified_maf0.01_popcornIN.tsv

score_file_eff=./precomputed_score/EUR_EAS_all_gen_eff.cscore
score_file_imp=/data1/jaeyoung/software/popcorn/EUR_EAS_all_gen_imp.cscore

fit_script=./popcorn/fit.py

popcorn fit -v 1 --cfile ${score_file_eff} \
                            --gen_effect \
                            --sfile1 ${s1} \
                            --sfile2 ${s2} \
                            ./popcorn_pge_${eur}_${eas}_maf0.01

popcorn fit -v 1 --cfile ${score_file_imp} \
                            --gen_effect \
                            --sfile1 ${s1} \
                            --sfile2 ${s2} \
                            ./popcorn_pgi_${eur}_${eas}_maf0.01

#!/bin/bash
source 00_global_variables.sh

MungeName=MetS_GWAS_munge
LDSCName=MetS_GWAS_ldsc

python2 ${src_LDSC} --h2 ${MungeName}.sumstats.gz \
                    --ref-ld-chr ./1000G_Phase3_ldscores/LDscore. \
                    --w-ld-chr ./1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                    --out "${LDSCName}"

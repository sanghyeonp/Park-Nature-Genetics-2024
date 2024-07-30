#!/bin/bash

SNPlist=w_hm3.snplist

SUMSTAT=MetS_GWAS.csv
EffN=1384348
MungeName=MetS_GWAS_munge

python2 ./ldsc/munge_sumstats.py --sumstats ${SUMSTAT} \
                        --N ${EffN} \
                        --out ${MungeName} \
                        --merge-alleles ${SNPlist} \
                        --chunksize 500000

#!/bin/bash
N_THREADS=8
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

Sumstat=QC.MetS_noUKB_GWAS.REFORMAT4PRS.GWAS.txt
SampleSize=1252787
OUT=./out;mkdir -p $OUT

ref_1kg_eur=ldblk_1kg_eur
UKB_bim=ukb_eur_unrel_comb

python3 ../PRScs/PRScs.py --ref_dir "${ref_1kg_eur}" \
                    --bim_prefix "${UKB_bim}" \
                    --sst_file "${Sumstat}" \
                    --n_gwas $SampleSize \
                    --phi 1e-2 \
                    --out_dir "${OUT}/MetS"

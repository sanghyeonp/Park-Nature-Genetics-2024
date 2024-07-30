#!/bin/bash
N_THREADS=8
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

Sumstat=MetS_noUKB_F1_GWAS.REFORMAT4PRS.GWAS.txt
SampleSize=279429
ref_1kg_eur=../ldblk_1kg_eur
UKB_bim=../ukb_eur_unrel_comb

OUT=./out_F1
mkdir -p ${OUT}

python3 ../PRScs/PRScs.py --ref_dir "${ref_1kg_eur}" \
                    --bim_prefix "${UKB_bim}" \
                    --sst_file "${Sumstat}" \
                    --n_gwas $SampleSize \
                    --phi 1e-2 \
                    --out_dir "${OUT}/MetS12_noUKB_F1"


Sumstat=MetS_noUKB_F2_GWAS.REFORMAT4PRS.GWAS.txt
SampleSize=286292
OUT=./out_F2
mkdir -p ${OUT}

python3 ../PRScs/PRScs.py --ref_dir "${ref_1kg_eur}" \
                    --bim_prefix "${UKB_bim}" \
                    --sst_file "${Sumstat}" \
                    --n_gwas $SampleSize \
                    --phi 1e-2 \
                    --out_dir "${OUT}/MetS12_noUKB_F2"


Sumstat=MetS_noUKB_F3_GWAS.REFORMAT4PRS.GWAS.txt
SampleSize=780473
OUT=./out_F3
mkdir -p ${OUT}

python3 ../PRScs/PRScs.py --ref_dir "${ref_1kg_eur}" \
                    --bim_prefix "${UKB_bim}" \
                    --sst_file "${Sumstat}" \
                    --n_gwas $SampleSize \
                    --phi 1e-2 \
                    --out_dir "${OUT}/MetS12_noUKB_F3"

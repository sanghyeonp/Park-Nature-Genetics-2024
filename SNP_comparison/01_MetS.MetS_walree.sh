#!/bin/bash
rscript=/snp_independence.R

trait1=MetS
trait2=MetS_walree

gwas1_snplist=../COJO/extract_cojo/cojo_out_jma.MetS.snplist.chr_pos
gwas2=../COJO/cojo_input/CF_MetS_results_adapt.txt.gz
r2_thres=0.01
p_thres=5e-8
window=500
n_thread=30

/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript ${rscript} \
    --gwas1-snplist ${gwas1_snplist} \
    --delim-gwas1-snplist whitespace \
    --snp-col-gwas1-snplist SNP \
    --chr-col-gwas1-snplist Chr \
    --pos-col-gwas1-snplist bp \
    --gwas2 ${gwas2} \
    --delim-gwas2 tab \
    --snp-col-gwas2 rsid \
    --chr-col-gwas2 chromosome \
    --pos-col-gwas2 base_pair_location \
    --pval-col-gwas2 p_value \
    --reference-panel "1kG" \
    --r2-threshold ${r2_thres} \
    --window ${window} \
    --pval-threshold ${p_thres} \
    --n-thread ${n_thread} \
    --prefix-out snp_independence.w_${window}.r2_${r2_thres}.p_${p_thres}.${trait1}_snplist.${trait2}_gwas \
    --delim-out comma

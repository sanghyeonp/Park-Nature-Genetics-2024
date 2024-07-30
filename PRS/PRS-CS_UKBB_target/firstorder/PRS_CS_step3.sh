#!/bin/bash
UKB_bim=ukb_eur_unrel_comb
plink --bfile ${UKB_bim} --score ./score_sum_F1.txt 2 4 6 sum --out ./score_sum_indiv_F1

plink --bfile ${UKB_bim} --score ./score_sum_F2.txt 2 4 6 sum --out ./score_sum_indiv_F2

plink --bfile ${UKB_bim} --score ./score_sum_F3.txt 2 4 6 sum --out ./score_sum_indiv_F3

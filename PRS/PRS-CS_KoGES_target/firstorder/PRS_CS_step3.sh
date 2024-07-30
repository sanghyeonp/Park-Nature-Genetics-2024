#!/bin/bash
KoGES_bim=KCHIP.snp.imputed.HRC.rsID.rmdup
plink --bfile ${KoGES_bim} --score ./score_sum_F1.txt 2 4 6 sum --out ./score_sum_indiv_F1

plink --bfile ${KoGES_bim} --score ./score_sum_F2.txt 2 4 6 sum --out ./score_sum_indiv_F2

plink --bfile ${KoGES_bim} --score ./score_sum_F3.txt 2 4 6 sum --out ./score_sum_indiv_F3

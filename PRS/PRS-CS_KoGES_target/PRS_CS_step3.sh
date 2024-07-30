#!/bin/bash

KoGES_bim=KCHIP.snp.imputed.HRC.rsID.rmdup
plink --bfile ${KoGES_bim} --score ./score_sum.txt 2 4 6 sum --out ./score_sum_indiv
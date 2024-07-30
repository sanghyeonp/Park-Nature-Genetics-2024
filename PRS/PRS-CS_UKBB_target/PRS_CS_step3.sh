#!/bin/bash

UKB_bim=ukb_eur_unrel_comb
plink --bfile ${UKB_bim} --score ./score_sum.txt 2 4 6 sum --out ./score_sum_indiv
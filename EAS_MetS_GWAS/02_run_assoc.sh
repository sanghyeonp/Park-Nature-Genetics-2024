#!/bin/bash

bfile=KCHIP.snp.imputed.HRC.abc
pheno_cov_file=KoGES_MetS_pheno_covar_table.txt
outputfile=KoGES_MetS_logistic_plink

# Phenotype name
pheno_name=MetS

# covariates
covnames=AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

###01.association test /w firth test
echo "start _________"
date

plink2  --bfile ${bfile} \
        --glm hide-covar firth-fallback \
        --ci 0.95 \
        --pheno ${pheno_cov_file} \
        --pheno-name ${pheno_name} \
        --covar ${pheno_cov_file} \
        --covar-name ${covnames} \
        --covar-variance-standardize \
        --allow-no-sex \
        --out ${outputfile} \
        --threads 16

echo "END___________"
date
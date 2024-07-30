library(data.table)
library(dplyr)


covariate_file <- "./UKB_covariates.csv"

df_cov <- fread(covariate_file, sep="\t", data.table = F, nThread = 1) %>%
    dplyr::select(f.eid, SEX, BIRTH_YEAR, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
    rename(FID = f.eid) %>%
    mutate(IID = FID) %>%
    dplyr::select(FID, IID, SEX, BIRTH_YEAR, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) 

df_cov[is.na(df_cov)] <- "NA"

write.table(df_cov,
            file = "ukb.cov.txt",
            sep=" ", row.names=F, quote=F)

mets_clinical_file <- "./UKB_MetS_clinical_phenotype.csv"

df_pheno <- fread(mets_clinical_file, sep=",", data.table = F, nThread = 1) %>%
    dplyr::select(f.eid, MetS_clinical) %>%
    rename(FID = f.eid,
            MetS = MetS_clinical) %>%
    mutate(IID = FID) %>%
    dplyr::select(FID, IID, MetS)

df_pheno[is.na(df_pheno)] <- "NA"

write.table(df_pheno,
            file = "ukb.pheno.txt",
            sep=" ", row.names=F, quote=F)
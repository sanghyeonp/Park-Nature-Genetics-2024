library(data.table)
library(dplyr)

### F1
trait <- "F1"
file_gwas <- "MetS_noUKB_F1_GWAS.csv"

df <- fread(file_gwas, sep=",", data.table=F, nThread = 20) %>%
    dplyr::select(SNP, CHR, BP, A1, A2, MAF, est, Pval_Estimate) %>%
    rename(BETA = est,
           P = Pval_Estimate)

write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)


### F2
trait <- "F2"
file_gwas <- "MetS_noUKB_F2_GWAS.csv"

df <- fread(file_gwas, sep=",", data.table=F, nThread = 20) %>%
    dplyr::select(SNP, CHR, BP, A1, A2, MAF, est, Pval_Estimate) %>%
    rename(BETA = est,
           P = Pval_Estimate)

write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### F3
trait <- "F3"
file_gwas <- "MetS_noUKB_F3_GWAS.csv"

df <- fread(file_gwas, sep=",", data.table=F, nThread = 20) %>%
    dplyr::select(SNP, CHR, BP, A1, A2, MAF, est, Pval_Estimate) %>%
    rename(BETA = est,
           P = Pval_Estimate)

write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)
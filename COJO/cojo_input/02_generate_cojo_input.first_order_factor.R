library(data.table)
library(dplyr)

df_bim <- fread("UKB_random10k_noqc.bim",
                sep="\t", data.table = F, col.names = c("chr", "snp", "cm", "bp", "ref_ukb", "alt_ukb"),
                nThread = 20) %>%
        dplyr::select(snp, ref_ukb)



### F1 GWAS
trait <- "F1"
file_gwas <- "MetS_F1_GWAS.withQsnp.csv"
df <- fread(file_gwas, sep=",", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, est, SE, Pval_Estimate) %>%
    mutate(N = 679472) %>%
    rename(freq = MAF,
        b = est,
        se = SE,
        p = Pval_Estimate)
        

df <- merge(df, df_bim, by.x = "SNP", by.y = "snp", all.x = T) %>%
    mutate(A1_new = ifelse(is.na(ref_ukb), A1,
                            ifelse(A1 == ref_ukb, A1, A2)),
            A2_new = ifelse(is.na(ref_ukb), A2,
                            ifelse(A1 == ref_ukb, A2, A1)),
            b_new = ifelse(is.na(ref_ukb), b,
                            ifelse(A1 == ref_ukb, b, -1 * b))
            ) %>%
    dplyr::select(SNP, A1_new, A2_new, freq, b_new, se, p, N) %>%
    rename(A1 = A1_new,
        A2 = A2_new,
        b = b_new)


write.table(df, 
            file = paste0("cojo_input.", trait, ".txt"),
            sep=" ", row.names=F, quote=F)

### F2 GWAS
trait <- "F2"
file_gwas <- "MetS_F2_GWAS.withQsnp.csv"
df <- fread(file_gwas, sep=",", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, est, SE, Pval_Estimate) %>%
    mutate(N = 728556) %>%
    rename(freq = MAF,
        b = est,
        se = SE,
        p = Pval_Estimate)
        

df <- merge(df, df_bim, by.x = "SNP", by.y = "snp", all.x = T) %>%
    mutate(A1_new = ifelse(is.na(ref_ukb), A1,
                            ifelse(A1 == ref_ukb, A1, A2)),
            A2_new = ifelse(is.na(ref_ukb), A2,
                            ifelse(A1 == ref_ukb, A2, A1)),
            b_new = ifelse(is.na(ref_ukb), b,
                            ifelse(A1 == ref_ukb, b, -1 * b))
            ) %>%
    dplyr::select(SNP, A1_new, A2_new, freq, b_new, se, p, N) %>%
    rename(A1 = A1_new,
        A2 = A2_new,
        b = b_new)


write.table(df, 
            file = paste0("cojo_input.", trait, ".txt"),
            sep=" ", row.names=F, quote=F)


### F3 GWAS
trait <- "F3"
file_gwas <- "MetS_F3_GWAS.withQsnp.csv"
df <- fread(file_gwas, sep=",", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, est, SE, Pval_Estimate) %>%
    mutate(N = 1086560) %>%
    rename(freq = MAF,
        b = est,
        se = SE,
        p = Pval_Estimate)
        

df <- merge(df, df_bim, by.x = "SNP", by.y = "snp", all.x = T) %>%
    mutate(A1_new = ifelse(is.na(ref_ukb), A1,
                            ifelse(A1 == ref_ukb, A1, A2)),
            A2_new = ifelse(is.na(ref_ukb), A2,
                            ifelse(A1 == ref_ukb, A2, A1)),
            b_new = ifelse(is.na(ref_ukb), b,
                            ifelse(A1 == ref_ukb, b, -1 * b))
            ) %>%
    dplyr::select(SNP, A1_new, A2_new, freq, b_new, se, p, N) %>%
    rename(A1 = A1_new,
        A2 = A2_new,
        b = b_new)


write.table(df, 
            file = paste0("cojo_input.", trait, ".txt"),
            sep=" ", row.names=F, quote=F)
library(data.table)
library(dplyr)

df_bim <- fread("UKB_random10k_noqc.bim",
                sep="\t", data.table = F, col.names = c("chr", "snp", "cm", "bp", "ref_ukb", "alt_ukb"),
                nThread = 20) %>%
        dplyr::select(snp, ref_ukb)



### MetS GWAS
file_gwas <- "MetS_GWAS.txt"
df <- fread(file_gwas, sep=" ", data.table=F, nThread=20) %>%
    dplyr::select(snp, a1, a2, freq, b, se, p, n) %>%
    rename(SNP = snp,
        A1 = a1,
        A2 = a2,
        N = n)

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
            file = "cojo_input.MetS.txt",
            sep=" ", row.names=F, quote=F)


### TG
file_gwas <- "CLEANED.TG_GLGC_UKB.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval, N) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval)

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
            file = "cojo_input.TG.txt",
            sep=" ", row.names=F, quote=F)

### HDL
file_gwas <- "CLEANED.HDL_GLGC_UKB.txt"

df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval, N) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval)

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
            file = "cojo_input.HDL.txt",
            sep=" ", row.names=F, quote=F)

### HTN
file_gwas <- "CLEANED.HTN_meta_FinngenR7_UKB.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval) %>%
    mutate(N = 508612)

write.table(df, 
            file = "cojo_input.HTN.txt",
            sep=" ", row.names=F, quote=F)

### FG
file_gwas <- "CLEANED.FG_MAGIC.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval, N) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval)

write.table(df, 
            file = "cojo_input.FG.txt",
            sep=" ", row.names=F, quote=F)

### T2D
file_gwas <- "CLEANED.T2D_meta_FinngenR7_Mahajan2022_MVP.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval) %>%
    mutate(N = 597437)

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
            file = "cojo_input.T2D.txt",
            sep=" ", row.names=F, quote=F)

### WC
file_gwas <- "CLEANED.WC_UKB.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval) %>%
    mutate(N = 385932)

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
            file = "cojo_input.WC.txt",
            sep=" ", row.names=F, quote=F)

### BMI
file_gwas <- "CLEANED.BMI_GIANT2018.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(SNP, A1, A2, MAF, Effect, SE, Pval, N) %>%
    rename(freq = MAF,
        b = Effect,
        se = SE,
        p = Pval)

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
            file = "cojo_input.BMI.txt",
            sep=" ", row.names=F, quote=F)

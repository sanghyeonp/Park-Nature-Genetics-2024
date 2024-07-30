library(data.table)
library(dplyr)

# IV list
file_in_iv <- "MetS_exposure_IV_list.csv"
df_iv <- fread(file_in_iv, sep=",", data.table = F, nThread = 10)
iv_list <- df_iv$SNP

# MetS no-UKB GWAS
file_in <- ".MetS_noUKB_GWAS.csv"

df <- fread(file_in, sep=",", data.table = F, nThread = 20) %>%
    filter(SNP %in% iv_list) %>%
    dplyr::select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate) %>%
    rename(POS = BP,
        maf_exp = MAF,
        a1_exp = A1,
        a2_exp = A2,
        b_exp = est,
        se_exp = SE,
        p_exp = Pval_Estimate)

write.table(df, 
            "exposure_data.MetS_noUKB.txt",
            sep=" ", row.names=F, col.names=T, quote=F)
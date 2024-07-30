library(data.table)
library(dplyr)

df_bim <- fread("EAS.bim",
                sep="\t", data.table = F, col.names = c("chr", "snp", "cm", "bp", "ref_1kg_eas", "alt_1kg_eas"),
                nThread = 20) %>%
        dplyr::select(snp, ref_1kg_eas)

### MetS GWAS KoGES
file_gwas <- "KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid.modified"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(RSID, effect_allele, other_allele, MAF, BETA, `LOG(OR)_SE`, P) %>%
    rename(SNP = RSID,
        A1 = effect_allele,
        A2 = other_allele,
        freq = MAF,
        b = BETA,
        se = `LOG(OR)_SE`,
        p = P) %>%
    mutate(N = 62314)

df <- merge(df, df_bim, by.x = "SNP", by.y = "snp", all.x = T) %>%
    mutate(A1_new = ifelse(is.na(ref_1kg_eas), A1,
                            ifelse(A1 == ref_1kg_eas, A1, A2)),
            A2_new = ifelse(is.na(ref_1kg_eas), A2,
                            ifelse(A1 == ref_1kg_eas, A2, A1)),
            b_new = ifelse(is.na(ref_1kg_eas), b,
                            ifelse(A1 == ref_1kg_eas, b, -1 * b))
            ) %>%
    dplyr::select(SNP, A1_new, A2_new, freq, b_new, se, p, N) %>%
    rename(A1 = A1_new,
        A2 = A2_new,
        b = b_new) %>%
    filter(nchar(SNP) != 0)

write.table(df, 
            file = "cojo_input.MetS.KoGES.txt",
            sep=" ", row.names=F, quote=F)
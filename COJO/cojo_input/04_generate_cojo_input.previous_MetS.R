library(data.table)
library(dplyr)


df_bim <- fread("UKB_random10k_noqc.bim",
                sep="\t", data.table = F, col.names = c("chr", "snp", "cm", "bp", "ref_ukb", "alt_ukb"),
                nThread = 20) %>%
        dplyr::select(snp, ref_ukb)



### MetS_lind GWAS
trait <- "MetS_lind"
file_gwas <- "UKBB_MetS_alla_Stefan.rsid.txt.gz"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(rsid, effect_allele, other_allele, eaf, beta, se, p_value, ntotal) %>%
    rename(SNP = rsid,
    A1 = effect_allele,
    A2 = other_allele,
    freq = eaf,
    b = beta,
    p = p_value,
    N = ntotal)

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

### MetS_walree GWAS
trait <- "MetS_walree"
file_gwas <- "CF_MetS_results_adapt.txt.gz"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=20) %>%
    dplyr::select(rsid, effect_allele, other_allele, minor_allele_frequency, beta, standard_error, p_value) %>%
    mutate(N = 461920) %>%
    rename(SNP = rsid,
    A1 = effect_allele,
    A2 = other_allele,
    freq = minor_allele_frequency,
    b = beta,
    se = standard_error,
    p = p_value)
        

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


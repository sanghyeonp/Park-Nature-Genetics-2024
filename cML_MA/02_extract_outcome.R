library(data.table)
library(dplyr)

# IV list
file_in_iv <- "MetS_exposure_IV_list.csv"
df_iv <- fread(file_in_iv, sep=",", data.table = F, nThread = 10)
iv_list <- df_iv$SNP

# Outcome list
file_in_mr_result <- "tsmr_best_sig_result.csv"
df_mr_result <- fread(file_in_mr_result, sep=",", data.table = F, nThread = 10,
                    colClasses = c("character", "character", "character",
                                "character", "integer", "numeric", "numeric", "numeric"))
mr_result.outcome_list <- df_mr_result$Phecode

# TSMR metadata
file_in_metadata <- "metadata_tsmr.tsv"
df_metadata <- fread(file_in_metadata, sep="\t", data.table = F, nThread = 10,
                    colClasses = c("character", "character", "character", "character", "character", "character")) %>%
                filter(Phecode %in% mr_result.outcome_list)


for (idx in 1:nrow(df_metadata)){
  phecode <- df_metadata[idx, ]$Phecode
  file_sumstat <- df_metadata[idx, ]$sumstat_path
  
  df_temp <- fread(file_sumstat, sep="\t", data.table = F, nThread = 20) %>%
    filter(rsid %in% iv_list) %>%
    dplyr::select(rsid, chr, pos, maf, alt, ref, beta, se, pval) %>%
    rename(SNP = rsid,
           CHR = chr,
           POS = pos,
           maf_out = maf,
           a1_out = alt,
           a2_out = ref,
           b_out = beta,
           se_out = se,
           p_out = pval)
  
  write.table(df_temp, 
              paste0("outcome_data.", phecode, ".txt"),
              sep=" ", row.names=F, col.names=T, quote=F)
}


library(data.table)
library(dplyr)
library(TwoSampleMR)

dir_out_har <- "tsmr_harmonized_data/"
dir.create(dir_out_har, showWarnings = FALSE)

dir_out_tsmr_res <- "tsmr_ivw_result/"
dir.create(dir_out_tsmr_res, showWarnings = FALSE)

# Exposure data
df_exp <- fread("exposure_data.MetS_noUKB.txt",
                sep=" ", data.table = F, nThread = 20)

df_tsmr_exp <- format_data(dat = df_exp,
                    type = "exposure",
                    snp_col = "SNP",
                    chr_col = "CHR",
                    pos_col = "POS",
                    beta_col = "b_exp",
                    se_col = "se_exp",
                    pval_col = "p_exp",
                    effect_allele_col = "a1_exp",
                    other_allele_col = "a2_exp"
                    )

# Metadata
file_in_mr_result <- "tsmr_best_sig_result.csv"
df_mr_result <- fread(file_in_mr_result, sep=",", data.table = F, nThread = 10,
                    colClasses = c("character", "character", "character",
                                "character", "integer", "numeric", "numeric", "numeric"))
mr_result.outcome_list <- df_mr_result$Phecode

# Outcome data
for (phecode in mr_result.outcome_list){
    file_in <- paste0("outcome_data.", phecode, ".txt")
    df_out <- fread(file_in, sep=" ", data.table = F, nThread = 20) 

    df_tsmr_out <- format_data(dat = df_out,
                    type = "outcome",
                    snp_col = "SNP",
                    chr_col = "CHR",
                    pos_col = "POS",
                    beta_col = "b_out",
                    se_col = "se_out",
                    pval_col = "p_out",
                    effect_allele_col = "a1_out",
                    other_allele_col = "a2_out"
                    )
    
    df_har <- harmonise_data(df_tsmr_exp, df_tsmr_out, action =2)
    
    saveRDS(df_har,
            file = paste0(dir_out_har, "/tsmr_harmonized.exp_MetSnoUKB.out_", phecode, ".rds"))
        
    tsmr_res <- mr(dat=df_har, method_list = c("mr_ivw"))
    saveRDS(tsmr_res,
            file = paste0(dir_out_tsmr_res, "/tsmr_ivw.exp_MetSnoUKB.out_", phecode, ".rds"))
}

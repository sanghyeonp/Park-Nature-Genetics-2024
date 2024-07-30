library(data.table)
library(dplyr)

dir_har <- "tsmr_harmonized_data/"

dir_out <- "cml_ma_input/"
dir.create(dir_out, showWarnings = FALSE)

# Metadata
file_in_mr_result <- "tsmr_best_sig_result.csv"
df_mr_result <- fread(file_in_mr_result, sep=",", data.table = F, nThread = 10,
                    colClasses = c("character", "character", "character",
                                "character", "integer", "numeric", "numeric", "numeric"))
mr_result.outcome_list <- df_mr_result$Phecode

# Outcome data
for (phecode in mr_result.outcome_list){
    df_har <- readRDS(paste0(dir_har, "tsmr_harmonized.exp_MetSnoUKB.out_", phecode, ".rds"))

    df_cml_input <- df_har %>%
                    dplyr::select(SNP, beta.exposure, se.exposure, beta.outcome, se.outcome) %>%
                    dplyr::rename(b_exp = beta.exposure,
                                se_exp = se.exposure,
                                b_out = beta.outcome,
                                se_out = se.outcome)


    saveRDS(df_cml_input,
            file = paste0(dir_out, "exp_MetSnoUKB.out_", phecode, ".rds"))
}


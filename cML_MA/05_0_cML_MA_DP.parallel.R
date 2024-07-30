sink("06_cML_MA_DP.parallel.log")

library(data.table)
library(dplyr)
library(MRcML)
library(foreach)
library(doParallel)

n_parallel <- 5

dir_in <- "cml_ma_input/"
dir_out <- "cml_ma_dp_result/"
dir.create(dir_out, showWarnings = FALSE)

# Metadata
file_in_mr_result <- "tsmr_best_sig_result.csv"
df_mr_result <- fread(file_in_mr_result, sep=",", data.table = F, nThread = 10,
                    colClasses = c("character", "character", "character",
                                "character", "integer", "numeric", "numeric", "numeric"))
mr_result.outcome_list <- df_mr_result$Phecode

# Define function
run_cML_MA_DP <- function(dir_in,
                        phecode,
                        dir_out
                        ){
    print(phecode)
    start_time <- Sys.time()
    df <- readRDS(paste0(dir_in, "exp_MetSnoUKB.out_", phecode, ".rds"))

    cML_dp_result <- mr_cML_DP(df$b_exp,
                    df$b_out,
                    df$se_exp,
                    df$se_out,
                    n = 1252786,
                    random_start = 20,
                    random_start_pert = 20,
                    random_seed = 12345,
                    num_pert = 200)
    
    saveRDS(cML_dp_result, 
            paste0(dir_out, "cML_dp_result.exp_MetSnoUKB.out_", phecode, ".rds"))

    end_time <- Sys.time()
    print(start_time)
    print(end_time)
    print(end_time - start_time)
}

# Register a parallel backend
cl <- makeCluster(n_parallel)
registerDoParallel(cl)

plink_ld_result <- foreach(phecode = mr_result.outcome_list,
                        .export = c("dir_in", "dir_out"),
                        .packages = c("MRcML")) %dopar% {
                    run_cML_MA_DP(dir_in = dir_in,
                        phecode = phecode,
                        dir_out = dir_out)
                    }
stopCluster(cl)
library(data.table)
library(dplyr)


dir.result <- "./cml_ma_dp_result"
result.list <- list.files(dir.result, pattern = "*.rds")

### Result to dataframe
df.combined <- NULL
for (idx in 1:length(result.list)){
    filename <- result.list[idx]
    phecode <- strsplit(gsub(".rds", "", filename), "_")[[1]][5]

    cML_dp_result <- readRDS(paste0(dir.result, "/", filename))
    df <- as.data.frame(t(data.frame(value = unlist(cML_dp_result))))
    rownames(df) <- NULL
    df$Phecode <- phecode
    
    if (idx == 1){
        df.combined <- df
    } else{
        df.combined <- df.combined %>%
            bind_rows(df)
    }
    
}

idx.list <- grep("*_invalid", colnames(df.combined))
df.bic_invalid <- df.combined[, idx.list]
invalid_bic <- c()
n_invalid_bic <- c()
for (idx in 1:nrow(df.bic_invalid)){
    row <- df.bic_invalid[idx, ] %>% 
        select_if(~ !any(is.na(.)))
    invalid_bic <- c(invalid_bic, paste(row, collapse = ";"))
    n_invalid_bic <- c(n_invalid_bic, length(colnames(row)))
    
}

df.combined2 <- df.combined[, -idx.list]
df.combined2$BIC_invalid <- invalid_bic
df.combined2$BIC_invalid_N <- n_invalid_bic

### Determine GOF
df.combined2 <- df.combined2 %>%
    mutate(approach_selection = ifelse(GOF1_p < 0.05 | GOF2_p < 0.05, "MA_BIC_DP", "MA_BIC"))


### Retain SNP that are invalid
dir_in <- "cml_ma_input/"
invalid_snp.list <- c()
for (idx in 1:nrow(df.combined2)){
    row <- df.combined2[idx,]
    phecode <- row$Phecode
    n_invalid <- row$BIC_invalid_N
    invalid_snp_idx <- row$BIC_invalid
    
    cML_input <- readRDS(paste0(dir_in, "exp_MetSnoUKB.out_", phecode, ".rds"))
    snp_list <- cML_input$SNP
    
    if (n_invalid == 0){
        invalid_snp <- ""
    } else{
        invalid_snp <- paste(snp_list[as.integer(strsplit(invalid_snp_idx, ";")[[1]])], collapse = ";")
    }
    
    invalid_snp.list <- c(invalid_snp.list, invalid_snp)
}

df.combined2$BIC_invalid_SNP <- invalid_snp.list


write.table(df.combined2,
            "table.06_1_combine_cML_MA_DP_result.csv",
            sep=",", row.names = F, quote = T)





library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

### 0. Data files
gwas_dis <- "MetS_GWAS.csv"
gwas_val <- "KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid.modified"

reference_panel <- "EUR"
reference_panel_chr <- "EUR_chr"
plink_exe <- "./plink/plink"
window <- 500
r2 <- 0.1
p1_thres <- 5e-8
p2_thres <- 5e-8

credible_window <- 50
credible_r2 <- 0.8

n_thread <- 20
### 1. Discovery GWAS에서 lead SNP 추리기.
df.gwas_dis <- fread(gwas_dis, sep=',', data.table = F, nThread = n_thread) %>%
    dplyr::select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate) %>%
    rename(POS = BP,
            a1_eur = A1,
           a2_eur = A2,
           af_eur = MAF,
            beta_eur = est,
           se_eur = SE,
           p_eur = Pval_Estimate)

file.plink_assoc <- "gwas_dis.assoc"
write.table(df.gwas_dis %>%
                dplyr::select(SNP, p_eur) %>%
                rename(P = p_eur), 
            file = file.plink_assoc, 
            sep = "\t", quote = F, row.names = F)

system(paste0(plink_exe, 
            " --bfile ", reference_panel, 
            " --clump ", file.plink_assoc, 
            " --clump-p1 ", p1_thres, 
            " --clump-p2 ", p2_thres, 
            " --clump-r2 ", r2, 
            " --clump-kb ", window, 
            " --out ", file.plink_assoc,
            " --threads ", n_thread), wait = TRUE)

df.leadsnp <- fread(paste0(file.plink_assoc, ".clumped"), data.table = F, nThread = n_thread) %>%
    dplyr::select(SNP) %>%
    rename(leadSNP_eur = SNP)

df.leadsnp <- merge(df.leadsnp, df.gwas_dis, by.x = "leadSNP_eur", by.y = "SNP", all.x = T)

### 2. credible set
compute_snp_ld <- function(snp1, snp2, plink_exe, file_reference, chr){
    plink_out <- system(paste0(plink_exe, 
                               " --bfile ", file_reference, chr, 
                               " --ld ", snp1, " ", snp2),
                        intern = TRUE,
                        wait = TRUE,
                        ignore.stderr = TRUE)
    ld_result <- grep("R-sq", plink_out, value = TRUE)
    
    if (identical(ld_result, character(0))){
        r2 <- "NA"
    } else{
        ld_result <- ld_result[1] # Refer to important note 1
        if (nchar(ld_result) != 0){
            ld_result_query <- regmatches(ld_result, regexpr("R-sq\\s+=\\s+([0-9.]+([eE][-+]?[0-9]+)?)", ld_result))
            
            r2 <- as.character(as.numeric(strsplit(ld_result_query, " ")[[1]][length(strsplit(ld_result_query, " ")[[1]])]))
            
        } else{
            r2 <- "NA"
        }
    }
    
    return (data.frame(ref_snp = snp1,
                       other_snp = snp2,
                       R2 = r2))
}


find_credible_set <- function(credible_idx,
                              df.leadsnp,
                              df_gwas,
                              leadsnp,
                              credible_window,
                              credible_r2,
                              plink_exe, reference_panel
){
    leadsnp <- df.leadsnp[credible_idx, ]$leadSNP_eur
    chr <- df_gwas[df_gwas$SNP == leadsnp, ]$CHR
    pos <- df_gwas[df_gwas$SNP == leadsnp, ]$POS
    pval <- df_gwas[df_gwas$SNP == leadsnp, ]$p_eur
    
    df_gwas.window <- df_gwas %>%
        filter(CHR == chr &
               POS >= pos - credible_window * 1000 &
               POS < pos + credible_window * 1000 &
               p_eur < pval * 100)
    
    if (nrow(df_gwas.window) > 0){
        snplist <- df_gwas.window$SNP
        df.r2 <- NULL
        for (snp in snplist){
            df.r2.single <- compute_snp_ld(snp1 = leadsnp, 
                                           snp2 = snp, 
                                           plink_exe = plink_exe, 
                                           file_reference = reference_panel, 
                                           chr = chr)
            df.r2 <- rbind(df.r2, df.r2.single)
        }
        
        df.r2 <- df.r2 %>%
            dplyr::select(-ref_snp) %>%
            filter(R2 != "NA") %>%
            mutate(R2 = as.numeric(R2))
        
        if (nrow(df.r2) > 0){
            df_gwas.window.r2 <- merge(df_gwas.window, df.r2, 
                                       by.x = "SNP", by.y = "other_snp", all.x = T) %>%
                filter(R2 >= credible_r2) %>%
                mutate(credible_set = credible_idx)
            return (df_gwas.window.r2)
        } else{
            return (data.frame())
        }
        
    } else{
        return (data.frame())
    }

}

# Register a parallel backend
cl <- makeCluster(n_thread)
registerDoParallel(cl)

df.credible <- foreach(credible_idx = 1:nrow(df.leadsnp),
                           .export = c("df.leadsnp", "df.gwas_dis", "credible_window", "credible_r2",
                                       "plink_exe", "reference_panel"),
                           .packages = c("dplyr", "data.table"),
                            .combine = "bind_rows") %dopar% {
                               
                           find_credible_set(credible_idx = credible_idx,
                                             df.leadsnp = df.leadsnp,
                                             df_gwas = df.gwas_dis,
                                             leadsnp = leadsnp,
                                             credible_window = credible_window,
                                             credible_r2 = credible_r2,
                                             plink_exe = plink_exe, 
                                             reference_panel = reference_panel_chr)
                           }

stopCluster(cl)

write.table(df.credible,
            "df_credible.csv",
            sep=",", row.names = F, quote = F)

df_credible <- fread("df_credible.csv", sep=",", data.table = F)


### 3. Validation GWAS information merge
df.gwas_val <- fread(gwas_val, sep="\t", data.table = F, nThread = n_thread) %>%
    dplyr::select(RSID, A1, MAF, BETA, `LOG(OR)_SE`, P) %>%
    rename(SNP = RSID,
           a1_eas = A1,
           freq_eas = MAF,
           beta_eas = BETA,
           se_eas = `LOG(OR)_SE`,
           p_eas = P)

df_credible2 <- merge(df_credible, df.gwas_val, by = "SNP", all.x = T) %>%
    mutate(a1_eas_new = ifelse(a1_eur == a1_eas, a1_eas, a1_eur),
           freq_eas_new = ifelse(a1_eur == a1_eas, freq_eas, 1 - freq_eas),
           beta_eas_new = ifelse(a1_eur == a1_eas, beta_eas, -1 * beta_eas)) %>%
    dplyr::select(-a1_eas, -freq_eas, -beta_eas) %>%
    rename(a1_eas = a1_eas_new,
           freq_eas = freq_eas_new,
           beta_eas = beta_eas_new) %>%
    na.omit()


df.credible.final <- NULL
for (credible_set_number in unique(df_credible2$credible_set)){
    print(credible_set_number)
    df_credible.filter <- df_credible2 %>% 
        filter(credible_set == credible_set_number &
                p_eas < 0.05 &
                ((beta_eur < 0 & beta_eas < 0) | (beta_eur > 0 & beta_eas > 0)))
    
    if (nrow(df_credible.filter) > 0){
        df_credible.filter <- df_credible.filter %>%
            arrange(p_eur)
        row_credible <- df_credible.filter[1, ] %>%
            mutate(is_leadsnp = SNP %in% df.leadsnp$leadSNP_eur)
    } else{
        row_credible <- data.frame()
    }
    df.credible.final <- rbind(df.credible.final, row_credible)
}


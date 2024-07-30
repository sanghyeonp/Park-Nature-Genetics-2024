library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(broom)
source("00_read_prs.R")


dir_out.coxph.tab <- "./coxph_PRS_table/"
dir.create(file.path(dir_out.coxph.tab), showWarnings = FALSE)

#########
file.df_merged <- "dataframe.coxph_survival_input.rds"
if (file.exists(file.df_merged)){
    df_merged <- readRDS(file.df_merged)
} else {
    ######### Read PRS
    df_prs <- read_PRS()

    ######### Read phenotype
    file_in.pheno_cvm <- "20240205_ukbb_377909_CVm_ID.txt"
    df_pheno.cvm <- fread(file_in.pheno_cvm, sep="\t", data.table = F, nThread = 5) %>%
    mutate(FID = wonlab_IID,
            IID = wonlab_IID) %>%
    dplyr::select(FID, IID, CVm, CVm_du_time, CVD_all, CVD_i, CVD_du_time)

    ######### Read covariates
    file_in.cov <- "reformat.ukb.covariates.csv"
    df_cov <- fread(file_in.cov, sep=",", data.table = F, nThread = 5)
    
    ######### Merge data
    # Merge phenotype and MetS PRS
    df_merged <- merge(df_pheno.cvm, df_prs, by = c("FID", "IID"), all.x = T)

    # Merge covariate
    df_merged <- merge(df_merged, df_cov, by = c("FID", "IID"), all.x = T)

    saveRDS(df_merged,
            file.df_merged)
}


################################# (Indcidence rate filter)
df_merged.incidence <- df_merged[!(df_merged$CVD_all == 1 & df_merged$CVD_i != 1), ]

surv_object1 <- Surv(time = df_merged.incidence$CVD_du_time, 
                    event = df_merged.incidence$CVD_i)

#################################
### MetS
prs <- "MetS"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.MetS + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

cox.zph(coxph_fit)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)

### BMI
prs <- "BMI"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.BMI + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)

### WC
prs <- "WC"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.WC + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)


### HTN
prs <- "HTN"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.HTN + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)
    

### FG
prs <- "FG"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.FG + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)

### T2D
prs <- "T2D"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.T2D + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)


### HDL
prs <- "HDL"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.HDL.rev + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)

### TG
prs <- "TG"
coxph_fit <- coxph(surv_object1 ~ PRS_scaled.TG + age_at_recruitment + factor(sex) +
                        genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                        genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                    data=df_merged.incidence)

write.table(as.data.frame(tidy(coxph_fit)) %>%
                mutate(PRS = prs,
                       estimate_lo_ci = estimate - 1.96 * std.error,
                       estimate_up_ci = estimate + 1.96 * std.error,
                       hr = exp(estimate),
                       hr_lo_ci = exp(estimate_lo_ci),
                       hr_up_ci = exp(estimate_up_ci)),
            paste0(dir_out.coxph.tab, "coxph_fit.PRS_", prs, ".csv"),
            sep=",", row.names = F, quote = F)
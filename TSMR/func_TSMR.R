library(TwoSampleMR)
library(MRInstruments)
library(MRPRESSO)
require(argparse)

func_run_MR <- function(exposure_name, outcome_name,
                    harmonized_rdsfile,
                    savedir_individual){
    
    dat0 <- readRDS(harmonized_rdsfile)

    no_MR <- FALSE

    if (is.null(dat0)){
        result <- rep("NA", each=16)
        result <- replace(result, c(1, 2, 3), c(outcome_name, exposure_name, "No IV left after harmonization."))
        no_MR <- TRUE
    }
    dat <- dat0[dat0$mr_keep == TRUE, ]
    if (nrow(dat) == 0){
        result <- rep("NA", each=16)
        result <- replace(result, c(1, 2, 3), c(outcome_name, exposure_name, "No IV left after harmonization and mr_keep."))
        no_MR <- TRUE
    }
    if (no_MR){
        colnames(result) <- c("Outcome", "Exposure", "Method", "N SNP", "Beta", "SE", "P-value", "Odds ratio", "OR L95", "OR U95",
                        "Heterogeneity Q", "Heterogeneity P-value", "MR-Egger intercept", "MR-Egger intercept P-value",
                        "MR-PRESSO global test P-value", "MR-PRESSO distortion test P-value")
        write.table(result, file = paste0(savedir_individual, "/", gsub(".harmonized.rds", ".csv", basename(harmonized_rdsfile))),
                sep=",", row.names=FALSE)
        break
    }

    # When only one IV is left
    if (nrow(dat) == 1){
        result <- rep("NA", each=16)
        result <- replace(result, c(1, 2, 3), c(outcome_name, exposure_name, "Only one IV present."))
        colnames(result) <- c("Outcome", "Exposure", "Method", "N SNP", "Beta", "SE", "P-value", "Odds ratio", "OR L95", "OR U95",
                        "Heterogeneity Q", "Heterogeneity P-value", "MR-Egger intercept", "MR-Egger intercept P-value",
                        "MR-PRESSO global test P-value", "MR-PRESSO distortion test P-value")
        write.table(result, file = paste0(savedir_individual, "/", gsub(".harmonized.rds", ".csv", basename(harmonized_rdsfile))),
                sep=",", row.names=FALSE)
        break
    }

    ## 4. MR
    # MR regression
    res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

    # Heterogeneity test
    heterogeneity_test <- mr_heterogeneity(dat, method_list = c("mr_ivw"))

    # Pleiotropy test: MR-Egger intercept for horizontal pleiotropy detection (p-value below 0.05 means there's outlier)
    pleiotropy_test <- mr_pleiotropy_test(dat)

    # MR-PRESSO
    ## MR-PRESSO 결과 해석: https://github.com/rondolab/MR-PRESSO/issues/13
    # res_mrpresso <- run_mr_presso(dat, NbDistribution = 1000) # Two-sample에서 제공해주는 MR-PRESSO 사용하지 말기 (Outlier가 있는데 MR-PRESSO 결과가 안 나올 때가 있음.)
    res_mrpresso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)

    res_final <- rbind(res, c(res$id.exposure[1], 
                            res$id.outcome[1],
                            res$outcome[1], 
                            res$exposure[1],
                            "MR-PRESSO", 
                            nrow(dat) - ifelse(is.na(res_mrpresso[[1]]$`Causal Estimate`[2]), 
                                                0, length(res_mrpresso[[2]]$`Distortion Test`$`Outliers Indices`)),
                            res_mrpresso[[1]]$`Causal Estimate`[2], 
                            res_mrpresso[[1]]$Sd[2],
                            res_mrpresso[[1]]$`P-value`[2]
                            ))

    res_final$outcome <- outcome_name
    res_final$exposure <- exposure_name

    res_final <- transform(res_final, id.exposure = as.character(id.exposure), 
                                    id.outcome = as.character(id.outcome),
                                    outcome = as.character(outcome),
                                    exposure = as.character(exposure),
                                    method = as.character(method),
                                    nsnp = as.numeric(nsnp),
                                    b = as.numeric(b),
                                    se = as.numeric(se),
                                    pval = as.numeric(pval))

    result <- generate_odds_ratios(res_final)

    result$Heterogeneity.Q <- heterogeneity_test$Q
    result$Heterogeneity.Pvalue <- heterogeneity_test$Q_pval
    result$MR_Egger.intercept <- pleiotropy_test$egger_intercept
    result$MR_Egger.Pvalue <- pleiotropy_test$pval
    if (is.na(res_mrpresso[[1]]$`Causal Estimate`[2])){
    result$MR_PRESSO.globaltest <- "NA"
    result$MR_PRESSO.distortiontest <- "NA"
    } else{
    result$MR_PRESSO.globaltest <- res_mrpresso[[2]]$`Global Test`$Pvalue
    result$MR_PRESSO.distortiontest <- res_mrpresso[[2]]$`Distortion Test`$Pvalue
    }
    

    result <- result[ , !(names(result) %in% c("id.exposure", "id.outcome", "lo_ci", "up_ci"))]

    colnames(result) <- c("Outcome", "Exposure", "Method", "N SNP", "Beta", "SE", "P-value", "Odds ratio", "OR L95", "OR U95",
                        "Heterogeneity Q", "Heterogeneity P-value", "MR-Egger intercept", "MR-Egger intercept P-value",
                        "MR-PRESSO global test P-value", "MR-PRESSO distortion test P-value")

    # saveRDS(result, file = paste0(savedir_individual, gsub(".MRinput.tsv", ".RDS", basename(exposure_file))))
    write.table(result, file = paste0(savedir_individual, "/", gsub(".harmonized.rds", ".csv", basename(harmonized_rdsfile))),
                sep=",", row.names=FALSE)


}


parser <- argparse::ArgumentParser(description=":: Run MR ::", formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--exposure_name", required=TRUE, 
                    help="")
parser$add_argument("--outcome_name", required=TRUE, 
                    help="")
parser$add_argument("--harmonized_rdsfile", required=TRUE, 
                    help="")
parser$add_argument("--savedir_individual", required=TRUE, 
                    help="")


## Get parser arguments
args <- parser$parse_args()
exposure_name <- args$exposure_name
outcome_name <- args$outcome_name
harmonized_rdsfile <- args$harmonized_rdsfile
savedir_individual <- args$savedir_individual


func_run_MR(exposure_name, outcome_name,
            harmonized_rdsfile,
            savedir_individual)

########
# func_run_MR(exposure_name="x", outcome_name="y",
                        # harmonized_rdsfile="/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/results/MR/MetS12_noUKB_MetS_exposure/harmonized_rds/981.gwas.imputed_v3.both_sexes.harmonized.rds",
                        # savedir_individual="/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/results/MR/MetS12_noUKB_MetS_exposure/TSMR_individual_result")
# ######## Global test significant
# exposure_name <- 'x'
# outcome_name <- 'y'
# harmonized_rdsfile <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/results/MR/MetS12_noUKB_MetS_exposure/harmonized_rds/20002_1223.gwas.imputed_v3.both_sexes.harmonized.rds"
# savedir_individual <- "NA"
# 
# dat0 <- readRDS(harmonized_rdsfile)
# dat <- dat0[dat0$mr_keep == TRUE, ]
# 
# ## 4. MR
# # MR regression
# res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
# 
# # Heterogeneity test
# heterogeneity_test <- mr_heterogeneity(dat, method_list = c("mr_ivw"))
# 
# # Pleiotropy test: MR-Egger intercept for horizontal pleiotropy detection (p-value below 0.05 means there's outlier)
# pleiotropy_test <- mr_pleiotropy_test(dat)
# 
# # MR-PRESSO
# res_mrpresso <- run_mr_presso(dat, NbDistribution = 400)
# res_mrpresso1 <- res_mrpresso
# 
# # Outlier test: MR-PRESSO global test (reject Ho, there is outlier) (p-value below 0.05 means there's outlier)
# outlier_test <- res_mrpresso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
# 
# res_final <- rbind(res, c(attributes(res_mrpresso)$id.exposure, 
#                           attributes(res_mrpresso)$id.outcome,
#                           attributes(res_mrpresso)$outcome, 
#                           attributes(res_mrpresso)$exposure,
#                           "MR-PRESSO", 
#                           nrow(dat) - length(res_mrpresso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
#                           res_mrpresso[[1]]$`Main MR results`$`Causal Estimate`, 
#                           res_mrpresso[[1]]$`Main MR results`$Sd,
#                           res_mrpresso[[1]]$`Main MR results`$`P-value`
# ))
# 
# res_final$outcome <- outcome_name
# res_final$exposure <- exposure_name
# 
# res_final <- transform(res_final, id.exposure = as.character(id.exposure), 
#                        id.outcome = as.character(id.outcome),
#                        outcome = as.character(outcome),
#                        exposure = as.character(exposure),
#                        method = as.character(method),
#                        nsnp = as.numeric(nsnp),
#                        b = as.numeric(b),
#                        se = as.numeric(se),
#                        pval = as.numeric(pval))
# 
# result <- generate_odds_ratios(res_final)
# 
# result$Heterogeneity.Q <- heterogeneity_test$Q
# result$Heterogeneity.Pvalue <- heterogeneity_test$Q_pval
# result$MR_Egger.intercept <- pleiotropy_test$egger_intercept
# result$MR_Egger.Pvalue <- pleiotropy_test$pval
# result$MR_PRESSO.globaltest <- outlier_test
# result$Causal <- if (sum(result$pval <= 0.05) >= 3) "Yes" else "No"
# 
# result$Best.causal.estimator <- best_method(result$Causal[1], 
#                                             result$Heterogeneity.Pvalue,
#                                             result$MR_Egger.Pvalue,
#                                             result$MR_PRESSO.globaltest
# )
# 
# best_method(result$Causal[1], 
#             result$Heterogeneity.Pvalue,
#             result$MR_Egger.Pvalue,
#             result$MR_PRESSO.globaltest
# )
# 
# result <- result[ , !(names(result) %in% c("id.exposure", "id.outcome", "lo_ci", "up_ci"))]
# 
# colnames(result) <- c("Outcome", "Exposure", "Method", "N SNP", "Beta", "SE", "P-value", "Odds ratio", "OR L95", "OR U95",
#                       "Heterogeneity Q", "Heterogeneity P-value", "MR-Egger intercept", "MR-Egger intercept P-value",
#                       "MR-PRESSO global test P-value", "Causal", "Best causal estimator")
# 
# result1 <- result
# 
# ######## Global test not significant
# exposure_name <- 'x1'
# outcome_name <- 'y2'
# harmonized_rdsfile <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/results/MR/MetS12_noUKB_MetS_exposure/harmonized_rds/M47.gwas.imputed_v3.both_sexes.harmonized.rds"
# savedir_individual <- "NA"
# 
# dat0 <- readRDS(harmonized_rdsfile)
# dat <- dat0[dat0$mr_keep == TRUE, ]
# 
# ## 4. MR
# # MR regression
# res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
# 
# # Heterogeneity test
# heterogeneity_test <- mr_heterogeneity(dat, method_list = c("mr_ivw"))
# 
# # Pleiotropy test: MR-Egger intercept for horizontal pleiotropy detection (p-value below 0.05 means there's outlier)
# pleiotropy_test <- mr_pleiotropy_test(dat)
# 
# # MR-PRESSO
# res_mrpresso <- run_mr_presso(dat, NbDistribution = 400)
# 
# # Outlier test: MR-PRESSO global test (reject Ho, there is outlier) (p-value below 0.05 means there's outlier)
# outlier_test2 <- res_mrpresso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
# 
# res_final <- rbind(res, c(attributes(res_mrpresso)$id.exposure, 
#                           attributes(res_mrpresso)$id.outcome,
#                           attributes(res_mrpresso)$outcome, 
#                           attributes(res_mrpresso)$exposure,
#                           "MR-PRESSO", 
#                           nrow(dat) - length(res_mrpresso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
#                           res_mrpresso[[1]]$`Main MR results`$`Causal Estimate`, 
#                           res_mrpresso[[1]]$`Main MR results`$Sd,
#                           res_mrpresso[[1]]$`Main MR results`$`P-value`
# ))
# 
# res_mrpresso[[1]]$`Main MR results`$`Causal Estimate`[2]
# 
# res_final$outcome <- outcome_name
# res_final$exposure <- exposure_name
# 
# res_final <- transform(res_final, id.exposure = as.character(id.exposure), 
#                        id.outcome = as.character(id.outcome),
#                        outcome = as.character(outcome),
#                        exposure = as.character(exposure),
#                        method = as.character(method),
#                        nsnp = as.numeric(nsnp),
#                        b = as.numeric(b),
#                        se = as.numeric(se),
#                        pval = as.numeric(pval))
# 
# result <- generate_odds_ratios(res_final)
# 
# result$Heterogeneity.Q <- heterogeneity_test$Q
# result$Heterogeneity.Pvalue <- heterogeneity_test$Q_pval
# result$MR_Egger.intercept <- pleiotropy_test$egger_intercept
# result$MR_Egger.Pvalue <- pleiotropy_test$pval
# result$MR_PRESSO.globaltest <- outlier_test
# result$Causal <- if (sum(result$pval <= 0.05) >= 3) "Yes" else "No"
# 
# result$Best.causal.estimator <- best_method(result$Causal[1], 
#                                             result$Heterogeneity.Pvalue,
#                                             result$MR_Egger.Pvalue,
#                                             result$MR_PRESSO.globaltest
# )
# 
# best_method(result$Causal[1], 
#             result$Heterogeneity.Pvalue,
#             result$MR_Egger.Pvalue,
#             result$MR_PRESSO.globaltest
# )
# 
# result <- result[ , !(names(result) %in% c("id.exposure", "id.outcome", "lo_ci", "up_ci"))]
# 
# colnames(result) <- c("Outcome", "Exposure", "Method", "N SNP", "Beta", "SE", "P-value", "Odds ratio", "OR L95", "OR U95",
#                       "Heterogeneity Q", "Heterogeneity P-value", "MR-Egger intercept", "MR-Egger intercept P-value",
#                       "MR-PRESSO global test P-value", "Causal", "Best causal estimator")

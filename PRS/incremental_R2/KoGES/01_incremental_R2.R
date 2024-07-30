library(data.table)
library(stringr)
library(boot)
library(DescTools)
set.seed(1)

prs_file <- c("../../PRS-CS_KoGES_target/score_sum_ukb_excluded_indiv.profile.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.BMI.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.WC.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.FG.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.T2D.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.HTN.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.HDL.reformat",
              "../../PRS-CS_KoGES_target/individual_phenotypes/score_sum_indiv.profile.TG.reformat",
              "../../PRS-CS_KoGES_target/firstorder/score_sum_indiv.profile.F1.reformat",
              "../../PRS-CS_KoGES_target/firstorder/score_sum_indiv.profile.F2.reformat",
              "../../PRS-CS_KoGES_target/firstorder/score_sum_indiv.profile.F3.reformat"
              )
prs_trait <- c("MetS", "BMI", "WC", "FG", "T2D", "HTN", "HDL", "TG", "F1", "F2", "F3")


idx <- 1
df_prs <- NULL
for (trait in prs_trait){
  file <- prs_file[idx]
  prs <- fread(file,
            sep="\t", nThread=2)
  prs2 <- prs[, c("IID", "SCORESUM")]
  colnames(prs2) <- c("IID", "PRS_raw")
  prs2$PRS_raw <- scale(prs2$PRS_raw)
  colnames(prs2) <- c("IID", paste0("PRS_", trait))
  if (idx == 1){
      df_prs <- prs2
  } else{
      df_prs <- merge(df_prs, prs2, by.x='IID', by.y='IID', all.x=TRUE, all.y=TRUE)
  }
  idx <- idx + 1
}

pheno <- fread("KoGES2022_MetS_binarization_phenotype_table.csv",
               sep=",", nThread=2)
pheno2 <- pheno[, c("DIST_ID", "Cohort", "AGE", "SEX", "MetS")]
colnames(pheno2) <- c("IID", "Cohort", "AGE", "SEX", "MetS")

data <- merge(pheno2, df_prs, by="IID")

pc <- fread("KCHIP.abc.snp.pca.eigenvec",
            sep="\t", nThread=2)
pc2 <- pc[, c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

data <- merge(data, pc2, by="IID")

print(head(data))


logistic_model <- function(target_phenotype, prs_trait, data, null_syntax, prs_syntax){
  nullmodel_syntax <- as.formula(null_syntax)
  prsmodel_syntax <- as.formula(prs_syntax)
  
  null_model <- glm(nullmodel_syntax, family=binomial(), data = data, na.action = na.omit)
  prs_model <- glm(prsmodel_syntax, family=binomial(), data = data, na.action = na.omit)
  
  prs_model.estimates <- summary(prs_model)$coefficients
  
  if (grepl("+", prs_trait, fixed = TRUE)){
    n_prs <- str_count(prs_trait,"\\+") + 1
    prs_estimate <- NULL
    prs_se <- NULL
    prs_pval <- NULL
    prs_estimate.upper_bound <- NULL
    prs_estimate.lower_bound <- NULL
    OR <- NULL
    OR_ci_lower <- NULL
    OR_ci_upper <- NULL
    for(i in 1:n_prs){
      prs_estimate_temp <- prs_model.estimates[i + 1, 1]
      prs_se_temp <- prs_model.estimates[i + 1, 2]
      prs_pval_temp <- prs_model.estimates[i + 1, 4]
      
      margin_of_error_temp <- prs_se_temp * qnorm(0.975)
      prs_estimate.upper_bound_temp <- prs_estimate_temp + margin_of_error_temp
      prs_estimate.lower_bound_temp <- prs_estimate_temp - margin_of_error_temp
      OR_temp <- exp(prs_estimate_temp)
      OR_ci_lower_temp <- exp(prs_estimate.lower_bound_temp)
      OR_ci_upper_temp <- exp(prs_estimate.upper_bound_temp)
      
      prs_estimate <- c(prs_estimate, prs_estimate_temp)
      prs_se <- c(prs_se, prs_se_temp)
      prs_pval <- c(prs_pval, prs_pval_temp)
      prs_estimate.upper_bound <- c(prs_estimate.upper_bound, prs_estimate.upper_bound_temp)
      prs_estimate.lower_bound <- c(prs_estimate.lower_bound, prs_estimate.lower_bound_temp)
      OR <- c(OR, OR_temp)
      OR_ci_lower <- c(OR_ci_lower, OR_ci_lower_temp)
      OR_ci_upper <- c(OR_ci_upper, OR_ci_upper_temp)
    }
    prs_estimate <- paste(prs_estimate, collapse=';')
    prs_se <- paste(prs_se, collapse=';')
    prs_pval <- paste(prs_pval, collapse=';')
    prs_estimate.upper_bound <- paste(prs_estimate.upper_bound, collapse=';')
    prs_estimate.lower_bound <- paste(prs_estimate.lower_bound, collapse=';')
    OR <- paste(OR, collapse=';')
    OR_ci_lower <- paste(OR_ci_lower, collapse=';')
    OR_ci_upper <- paste(OR_ci_upper, collapse=';')
  }else{
    prs_estimate <- prs_model.estimates[2, 1]
    prs_se <- prs_model.estimates[2, 2]
    prs_pval <- prs_model.estimates[2, 4]
    
    margin_of_error <- prs_se * qnorm(0.975)
    prs_estimate.upper_bound <- prs_estimate + margin_of_error
    prs_estimate.lower_bound <- prs_estimate - margin_of_error
    
    OR <- exp(prs_estimate)
    OR_ci_lower <- exp(prs_estimate.lower_bound)
    OR_ci_upper <- exp(prs_estimate.upper_bound)
  }
  
  null.r2 <- PseudoR2(null_model, which = r2_eval)
  prs.r2 <- PseudoR2(prs_model, which = r2_eval)
  incremental.r2 <- prs.r2 - null.r2
  
  n_observation <- nobs(prs_model)
  return (cbind(target_phenotype, prs_trait, n_observation, prs_estimate, prs_se, prs_pval, prs_estimate.lower_bound, prs_estimate.upper_bound, 
                OR, OR_ci_lower, OR_ci_upper, 
                null.r2, prs.r2, incremental.r2))
}

bs_fnc <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  
  null_formula <- as.formula(str_split(formula, ";")[[1]][1])
  null.fit <- glm(null_formula, data=d, family=binomial(), na.action = na.omit)
  r2.null <- PseudoR2(null.fit, which = r2_eval)
  
  prs_formula <- as.formula(str_split(formula, ";")[[1]][2])
  prs.fit <- glm(prs_formula, data=d, family=binomial(), na.action = na.omit)
  prs_model.estimates <- summary(prs.fit)$coefficients
  estimate <- prs_model.estimates[2, 1]
  r2.prs <- PseudoR2(prs.fit, which = r2_eval)
  
  incremental_r2 <- r2.prs - r2.null
  
  return(incremental_r2)
}

r2_eval <- "McFaddenAdj"

n_boot <- 1000

bs_ci_type <- "perc"
target_phenotype <- "MetS"


results <- NULL

for (trait in prs_trait){
  prs <- paste0("PRS_", trait)
  cat(paste0("Running --> PRS: ", prs, " / Target phenotype: ", target_phenotype, "\n"))

  null_syntax <- paste0(target_phenotype, " ~", " SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  prs_syntax <- paste0(target_phenotype, " ~", " ", prs, " + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  full_syntax <- paste(null_syntax, prs_syntax, sep=";")

  result <- logistic_model(target_phenotype, prs, data, null_syntax, prs_syntax)
  # Bootstrap to obtain CI for the incremental R2
  boot_result <- boot(data=data, statistic=bs_fnc, R=n_boot, formula=full_syntax)
  boot_ci <- boot.ci(boot_result, type=bs_ci_type, index=1)

  incrementalR2.lower_bound <- boot_ci$percent[length(boot_ci$percent) - 1]
  incrementalR2.upper_bound <- boot_ci$percent[length(boot_ci$percent)]

  result <- cbind(result, incrementalR2.lower_bound, incrementalR2.upper_bound)
  cat(result)
  cat('\n')
  results <- rbind(results, result)
}


results <- as.data.frame(results)
print(results)

results <- setNames(results, c("Trait", "PRS", "N_obs", "PRS estimate", "SE", "P", "PRS estimate lower bound", "PRS estimate upper bound", 
                              "OR", "OR lower bound", "OR upper bound", 
                              "Null R2", "PRS R2", "Incremental R2 (%)", "Incremental R2 (%) lower bound", "Incremental R2 (%) upper bound"))

write.csv(results, file="KoGES_incremental_R2.csv", row.names=FALSE)
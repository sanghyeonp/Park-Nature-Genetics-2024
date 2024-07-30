
source("./01_load_data.R")
library(stringr)
library(boot)
library(DescTools)
set.seed(1)

save.dir <- "./"
result_file_name <- "UKB_incremental_R2.csv"

PRS.names <- c("PRS_BMI", "PRS_FG", "PRS_HDL", "PRS_HTN", "PRS_T2D", 
                "PRS_TG", "PRS_WC", "PRS_F1", "PRS_F2", "PRS_F3", "PRS_MetS")


phenotype.names <- c("MetS_clinical")

r2_eval <- "McFaddenAdj"

n_boot <- 1000
bs_ci_type <- "perc"

###########################################################################
data <- load_all_data()
data <- data[!(is.na(data$MetS_clinical)) & !(is.na(data$PRS_MetS)), ]
data$MetS_clinical <- as.factor(data$MetS_clinical)
data$SEX <- as.factor(data$SEX)



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
  d <- data[indices,]
  
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


results <- NULL

for(prs in PRS.names){
  for(target_phenotype in phenotype.names){
    cat(paste0("Running --> PRS: ", prs, " / Target phenotype: ", target_phenotype, "\n"))

    null_syntax <- paste0(target_phenotype, " ~", " SEX + BIRTH_YEAR + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    prs_syntax <- paste0(target_phenotype, " ~", " ", prs, " + SEX + BIRTH_YEAR + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    full_syntax <- paste(null_syntax, prs_syntax, sep=";")

    if(target_phenotype %in% quantitative){
      do_later <- NULL
    }else{
      result <- logistic_model(target_phenotype, prs, data, null_syntax, prs_syntax)
      # Bootstrap to obtain CI for the incremental R2
      boot_result <- boot(data=data, statistic=bs_fnc, R=n_boot, formula=full_syntax)
      boot_ci <- boot.ci(boot_result, type=bs_ci_type, index=1)

      incrementalR2.lower_bound <- boot_ci$percent[length(boot_ci$percent) - 1]
      incrementalR2.upper_bound <- boot_ci$percent[length(boot_ci$percent)]

      result <- cbind(result, incrementalR2.lower_bound, incrementalR2.upper_bound)
      cat(result)
      cat('\n')
    }
    results <- rbind(results, result)
  }
}

results <- as.data.frame(results)
print(results)

results <- setNames(results, c("Trait", "PRS", "N_obs", "PRS estimate", "SE", "P", "PRS estimate lower bound", "PRS estimate upper bound", 
                                "OR", "OR lower bound", "OR upper bound", 
                                "Null R2", "PRS R2", "Incremental R2 (%)", "Incremental R2 (%) lower bound", "Incremental R2 (%) upper bound"))

write.csv(results, file=result_file_name, row.names=FALSE)

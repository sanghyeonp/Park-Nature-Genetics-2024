
cat("Reading PRS...\n")
prs_file <- "../../PRS-CS/firstorder/PRS_MetS_F3.csv"
df1 <- read.csv(prs_file, sep=',', header=TRUE)
colnames(df1) <- c('FID', 'eid', 'PHENO', 'CNT', 'CNT2', 'PRS')
df1 <- df1[, c("eid", "PRS")]
df1$PRS <- scale(df1$PRS) # Normalize PRS

cat("Reading covariates...\n")
covariates_file <- "../UKB_covariates.csv"
df2 <- read.csv(covariates_file, header = TRUE, sep='\t')

cat("Merging PRS and covariates...\n")
df <- merge(df1, df2, by.x='eid', by.y='f.eid', all.x=TRUE, all.y=TRUE)

cat("Reading Phecodes...\n")
phecode_file <- "../UKB_Phecode_new"
df3 <- read.csv(phecode_file, header = TRUE, sep='\t')

cat("Merging Phecodes...\n")
data <- merge(df, df3, by.x='eid', by.y='eid', all.x=TRUE, all.y=TRUE)


phenotypes <- colnames(data)[grepl("phecode_", colnames(data))]

n_phecode <- length(phenotypes)

z_95 <- qnorm(0.975,mean=0,sd=1)
idx <- 1
results <- NULL
for(phecode in phenotypes){
  prsmodel_syntax <- as.formula(paste0(phecode, " ~", " PRS + SEX + BIRTH_YEAR + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  prs_model <- glm(prsmodel_syntax, family=binomial(), data = data, na.action = na.omit)
  
  n_observation <- nobs(prs_model)
  beta <- as.numeric(coef(summary(prs_model))[,'Estimate'][2])
  SE <- as.numeric(coef(summary(prs_model))[,'Std. Error'][2])
  beta_ci_lower <- beta - (z_95 * SE)
  beta_ci_upper <- beta + (z_95 * SE)
  Z <- as.numeric(coef(summary(prs_model))[,'z value'][2])
  P <- as.numeric(coef(summary(prs_model))[,'Pr(>|z|)'][2])
  OR <- exp(beta)
  OR_ci_lower <- exp(beta_ci_lower)
  OR_ci_upper <- exp(beta_ci_upper)
    
  result <- cbind(phecode, n_observation, beta, SE, beta_ci_lower, beta_ci_upper, Z, P, OR, OR_ci_lower, OR_ci_upper)
  
  idx <- idx + 1
  cat(paste0("[", idx, "/", n_phecode, "]\n"))
  cat(result)
  cat('\n')
  results <- rbind(results, result)
}


results <- as.data.frame(results)

results <- setNames(results, c("Phecode", "n_obs", "Beta", "SE", "Beta_CI_lower", "Beta_CI_upper", "Z-score", "P-value", 
                               "OR", "OR_CI_lower", "OR_CI_upper"))

write.csv(results, file="MetS_noUKB_F3_PRSPheWAS_results.csv")

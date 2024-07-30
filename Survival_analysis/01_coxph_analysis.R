library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
source("00_read_prs.R")

######### References

######### Read PRS
df_prs <- read_PRS()

######### Read covariates
file_in.cov <- "reformat.ukb.covariates.csv"
df_cov <- fread(file_in.cov, sep=",", data.table = F, nThread = 5)

######### Read phenotype
file_in.pheno_cvm <- "20240205_ukbb_377909_CVm_ID.txt"
df_pheno.cvm <- fread(file_in.pheno_cvm, sep="\t", data.table = F, nThread = 5) %>%
       dplyr::select(-CAD_PRS) %>%
       dplyr::select(FID, IID, CVm, CVm_du_time, CVD_all, CVD_i, CVD_du_time)


######### Merge data
# Merge phenotype and MetS PRS
df_merged <- merge(df_pheno.cvm, df_prs.mets, by = c("FID", "IID"), all.x = T)

# Merge covariate
df_merged <- merge(df_merged, df_cov, by = c("FID", "IID"), all.x = T)

######### Cox regression (univariate)
### uni_model1: MetS PRS (all)
res.cox.uni_model1 <- coxph(Surv(CV_time_month, CV) ~ PRS_scaled.MetS, data = df_merged)
summary(res.cox.uni_model1)

plot.res.cox.uni_model1 <- ggsurvplot(survfit(res.cox.uni_model1), 
           data = df_merged,
           conf.int = TRUE,
           palette = "#2E9FDF",
           xlab = "Months",
           ggtheme = theme_minimal())


pdf("res.cox.uni_model1.pdf")
print(plot.res.cox.uni_model1, 
      newpage = FALSE)
dev.off()
ggsave("res.cox.uni_model1.pdf",
       plot.res.cox.uni_model1, device = "pdf", scale = 2,
       units = "mm", width = 50, height = 50)

### uni_model2: MetS PRS (categorized)
res.cox.uni_model2 <- coxph(Surv(CV_time_month, CV) ~ factor(PRS_group.MetS), data = df_merged)
summary(res.cox.uni_model2)

######### Cox regression (multivariate)
### 1. multi_model1: MetS PRS, age, sex, genetic PC 1-10
res.cox.multi_model1 <- coxph(Surv(CV_time, CV) ~ PRS_scaled.MetS + age_at_recruitment + factor(sex) +
                                  genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                                  genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                              data = df_merged)
summary(res.cox.multi_model1)

### 2. multi_model2: categorized MetS PRS, age, sex, genetic PC 1-10
res.cox.multi_model2 <- coxph(Surv(CV_time, CV) ~ factor(PRS_group.MetS) + age_at_recruitment + factor(sex) +
                                  genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                                  genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10, 
                              data = df_merged)
summary(res.cox.multi_model2)

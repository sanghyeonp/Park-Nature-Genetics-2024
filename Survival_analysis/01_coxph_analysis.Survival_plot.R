
library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(broom)
source("../00_read_prs.R")


######## Pre-defined functions
pretty_number <- function(x, n_digit){
        return(prettyNum(round(x, n_digit), big.mark=",", scientific=FALSE, drop0trailing=FALSE))
}

######## Make directory
dir_out.coxph_group.rds <- "./coxph_PRS_group_rds/"
dir.create(file.path(dir_out.coxph_group.rds), showWarnings = FALSE)

dir_out.coxph_group.tab <- "./coxph_PRS_group_table/"
dir.create(file.path(dir_out.coxph_group.tab), showWarnings = FALSE)

dir_out.survfit_plot <- "./survfit_plot/"
dir.create(file.path(dir_out.survfit_plot), showWarnings = FALSE)


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

coxph_fit.mets <- coxph(surv_object1 ~ factor(PRS_group.MetS) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.mets,
        rds.coxph_fit)

df.coxph_fit.mets <- as.data.frame(tidy(coxph_fit.mets)) %>%
    mutate(PRS = prs,
           estimate_lo_ci = estimate - 1.96 * std.error,
           estimate_up_ci = estimate + 1.96 * std.error,
           hr = exp(estimate),
           hr_lo_ci = exp(estimate_lo_ci),
           hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.mets,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.mets  <- survfit(surv_object1 ~ factor(PRS_group.MetS), 
                                data=df_merged.incidence)

ggsurv_plot.mets <- ggsurvplot(surv_fit.mets, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.mets$plot <- ggsurv_plot.mets$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)2"))$hr, 4),
                            " (", 
                            pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.mets %>% filter(term == "factor(PRS_group.MetS)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.mets, newpage = FALSE)
dev.off()

### BMI
prs <- "BMI"

coxph_fit.bmi <- coxph(surv_object1 ~ factor(PRS_group.BMI) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.bmi,
        rds.coxph_fit)

df.coxph_fit.bmi <- as.data.frame(tidy(coxph_fit.bmi)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.bmi,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.bmi  <- survfit(surv_object1 ~ factor(PRS_group.BMI), 
                                data=df_merged.incidence)

ggsurv_plot.bmi <- ggsurvplot(surv_fit.bmi, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.bmi$plot <- ggsurv_plot.bmi$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.bmi %>% filter(term == "factor(PRS_group.BMI)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.bmi, newpage = FALSE)
dev.off()

### WC
prs <- "WC"
coxph_fit.wc <- coxph(surv_object1 ~ factor(PRS_group.WC) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.wc,
        rds.coxph_fit)

df.coxph_fit.wc <- as.data.frame(tidy(coxph_fit.wc)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.wc,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.wc  <- survfit(surv_object1 ~ factor(PRS_group.WC), 
                                data=df_merged.incidence)

ggsurv_plot.wc <- ggsurvplot(surv_fit.wc, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.wc$plot <- ggsurv_plot.wc$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.wc %>% filter(term == "factor(PRS_group.WC)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.wc, newpage = FALSE)
dev.off()


### HTN
prs <- "HTN"
coxph_fit.htn <- coxph(surv_object1 ~ factor(PRS_group.HTN) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.htn,
        rds.coxph_fit)

df.coxph_fit.htn <- as.data.frame(tidy(coxph_fit.htn)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.htn,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.htn  <- survfit(surv_object1 ~ factor(PRS_group.HTN), 
                                data=df_merged.incidence)

ggsurv_plot.htn <- ggsurvplot(surv_fit.htn, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.htn$plot <- ggsurv_plot.htn$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.htn %>% filter(term == "factor(PRS_group.HTN)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.htn, newpage = FALSE)
dev.off()


### T2D
prs <- "T2D"
coxph_fit.t2d <- coxph(surv_object1 ~ factor(PRS_group.T2D) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.t2d,
        rds.coxph_fit)

df.coxph_fit.t2d <- as.data.frame(tidy(coxph_fit.t2d)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.t2d,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.t2d  <- survfit(surv_object1 ~ factor(PRS_group.T2D), 
                                data=df_merged.incidence)

ggsurv_plot.t2d <- ggsurvplot(surv_fit.t2d, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.t2d$plot <- ggsurv_plot.t2d$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.t2d %>% filter(term == "factor(PRS_group.T2D)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.t2d, newpage = FALSE)
dev.off()


### FG
prs <- "FG"
coxph_fit.fg <- coxph(surv_object1 ~ factor(PRS_group.FG) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.fg,
        rds.coxph_fit)

df.coxph_fit.fg <- as.data.frame(tidy(coxph_fit.fg)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.fg,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.fg  <- survfit(surv_object1 ~ factor(PRS_group.FG), 
                                data=df_merged.incidence)

ggsurv_plot.fg <- ggsurvplot(surv_fit.fg, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.fg$plot <- ggsurv_plot.fg$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.fg %>% filter(term == "factor(PRS_group.FG)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.fg, newpage = FALSE)
dev.off()

### HDL
prs <- "HDL"
coxph_fit.hdl <- coxph(surv_object1 ~ factor(PRS_group.HDL) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.hdl,
        rds.coxph_fit)

df.coxph_fit.hdl <- as.data.frame(tidy(coxph_fit.hdl)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.hdl,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.hdl  <- survfit(surv_object1 ~ factor(PRS_group.HDL), 
                                data=df_merged.incidence)

ggsurv_plot.hdl <- ggsurvplot(surv_fit.hdl, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.hdl$plot <- ggsurv_plot.hdl$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.hdl %>% filter(term == "factor(PRS_group.HDL)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.hdl, newpage = FALSE)
dev.off()


### HDL - Note that PRS of HDL has been reversed
prs <- "HDL.rev"
coxph_fit.hdl.rev <- coxph(surv_object1 ~ factor(PRS_group.HDL.rev) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.hdl.rev,
        rds.coxph_fit)

df.coxph_fit.hdl.rev <- as.data.frame(tidy(coxph_fit.hdl.rev)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.hdl.rev,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.hdl.rev  <- survfit(surv_object1 ~ factor(PRS_group.HDL.rev), 
                                data=df_merged.incidence)

ggsurv_plot.hdl.rev <- ggsurvplot(surv_fit.hdl.rev, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.hdl.rev$plot <- ggsurv_plot.hdl.rev$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.hdl.rev %>% filter(term == "factor(PRS_group.HDL.rev)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.hdl.rev, newpage = FALSE)
dev.off()

### TG
prs <- "TG"
coxph_fit.tg <- coxph(surv_object1 ~ factor(PRS_group.TG) + age_at_recruitment + factor(sex) +
                            genetic_PC1 + genetic_PC2 + genetic_PC3 + genetic_PC4 + genetic_PC5 +
                            genetic_PC6 + genetic_PC7 + genetic_PC8 + genetic_PC9 + genetic_PC10,
                        data = df_merged.incidence)

rds.coxph_fit <- paste0(dir_out.coxph_group.rds, "coxph_fit.PRS_", prs, ".4groups.rds")
saveRDS(coxph_fit.tg,
        rds.coxph_fit)

df.coxph_fit.tg <- as.data.frame(tidy(coxph_fit.tg)) %>%
    mutate(PRS = prs,
        estimate_lo_ci = estimate - 1.96 * std.error,
        estimate_up_ci = estimate + 1.96 * std.error,
        hr = exp(estimate),
        hr_lo_ci = exp(estimate_lo_ci),
        hr_up_ci = exp(estimate_up_ci))

write.table(df.coxph_fit.tg,
            paste0(dir_out.coxph_group.tab, "coxph_fit.PRS_", prs, ".4groups.csv"),
            sep=",", row.names = F, quote = F)


surv_fit.tg  <- survfit(surv_object1 ~ factor(PRS_group.TG), 
                                data=df_merged.incidence)

ggsurv_plot.tg <- ggsurvplot(surv_fit.tg, 
                        # plot cumulative events
                        fun = "event", 
                        # plot start from origin
                        axes.offset = TRUE,
                        legend.labs=c("Low-risk", 
                                    "Intermediate-risk",
                                    "High-risk",
                                    "Very high-risk"),
                        palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
                        # break by 2 years
                        xlim = c(0, 10), break.time.by = 2,  
                        ylim = c(0, 0.16), 
                        # No median survival
                        surv.median.line="none",
                        # Show confidence interval
                        conf.int = F, 
                        pval=F, 
                        # Number at risk and cumulative number of evenets
                        risk.table="nrisk_cumevents", 
                        # Do not draw censors
                        censor = FALSE,
                        # Y-axis as percentage
                        surv.scale = c("percent"),
                        xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
                        legend.title="", legend = c("top"),
                        font.x = c(17, "bold", "black"),
                        font.y = c(17, "bold", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        fontsize = 4,
                        size = 1, 
                        ggtheme = theme_classic())

ggsurv_plot.tg$plot <- ggsurv_plot.tg$plot +
    ggplot2::annotate("text", 
                    x = 0, y = 0.15, hjust = 0, vjust = 1,# x and y coordinates of the text
                    label = paste0("Low-risk; Reference\n", 
                    "Intermediate-risk; HR = ", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)2"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)2"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)2"))$hr_up_ci, 4),
                            ")\n", 
                    "High-risk; HR = ", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)3"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)3"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)3"))$hr_up_ci, 4),
                            ")\n", 
                    "Very high-risk; HR = ", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)4"))$hr, 4),
                            " (", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)4"))$hr_lo_ci, 4),
                            "-", pretty_number((df.coxph_fit.tg %>% filter(term == "factor(PRS_group.TG)4"))$hr_up_ci, 4),
                            ")\n"), 
                    size = 5) +
    theme(legend.text = element_text(size = 14, color = "black"))


pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_", prs, "_4groups.pdf"), 
    width = 13, height = 11)
print(ggsurv_plot.tg, newpage = FALSE)
dev.off()

# ###############################
# splots <- list()
# splots[[1]] <- ggsurv_plot.mets
# splots[[2]] <- ggsurv_plot.bmi
# splots[[3]] <- ggsurv_plot.wc
# splots[[4]] <- ggsurv_plot.htn
# splots[[5]] <- ggsurv_plot.fg
# splots[[6]] <- ggsurv_plot.t2d
# splots[[7]] <- ggsurv_plot.hdl
# splots[[8]] <- ggsurv_plot.tg
# 
# pdf(paste0(dir_out.survfit_plot, "CVD_incidence_rate.PRS_combined_4groups.pdf"), 
#     width = 18, height = 22.5)
# arrange_ggsurvplots(splots, print = TRUE,
#                     ncol = 2, nrow = 4)
# dev.off()
# 
# ###############################
# require(devtools)
# devtools::install_github("adamleejohnson/R-ajtools")
# library(ajtools)
# ggsurv_plot.tg2 <- ajtools::ggsurvplot_patchwork(surv_fit.tg, 
#                              # plot cumulative events
#                              fun = "event", 
#                              # plot start from origin
#                              axes.offset = TRUE,
#                              legend.labs=c("Low-risk", 
#                                            "Intermediate-risk",
#                                            "High-risk",
#                                            "Very high-risk"),
#                              palette = c("#002663", "#009f4d", "#efdf00", "#e4002b"),
#                              # break by 2 years
#                              xlim = c(0, 10), break.time.by = 2,  
#                              ylim = c(0, 0.16), 
#                              # No median survival
#                              surv.median.line="none",
#                              # Show confidence interval
#                              conf.int = F, 
#                              pval=F, 
#                              # Number at risk and cumulative number of evenets
#                              risk.table="nrisk_cumevents", 
#                              # Do not draw censors
#                              censor = FALSE,
#                              # Y-axis as percentage
#                              surv.scale = c("percent"),
#                              xlab = "Years", ylab = "Cumulative CVD incidence rate (%)",
#                              legend.title="", legend = c("top"),
#                              font.x = c(14, "bold", "black"),
#                              font.y = c(14, "bold", "black"),
#                              font.tickslab = c(12, "plain", "black"),
#                              fontsize = 4,
#                              size = 1, 
#                              ggtheme = theme_classic())
# 
# ggsurv_plot.tg2 + ggsurv_plot.tg2 + plot_layout(nrow = 4, byrow = FALSE)

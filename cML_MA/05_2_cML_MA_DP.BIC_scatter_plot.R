library(data.table)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(scales)

dir_cml_input <- "cml_ma_input/"
dir_res <- "cml_ma_result/"
dir_tsmr_res <- "tsmr_ivw_result/"
dir_out <- "cml_ma_dp_result.cML_IVW_scatter_plot/"
dir.create(dir_out, showWarnings = FALSE)

# cML results
file_in_cml <- "table.06_1_combine_cML_MA_DP_result.csv"
df.cml <- fread(file_in_cml, sep=",", data.table = F, nThread = 1,
                    colClasses = c(rep("numeric", 14), rep("character", 2), "integer", rep("character", 2)))

# Metadata
file_in_mr_result <- "tsmr_best_sig_result.csv"
df_mr_result <- fread(file_in_mr_result, sep=",", data.table = F, nThread = 10,
                    colClasses = c("character", "character", "character",
                                "character", "integer", "numeric", "numeric", "numeric"))


# TSMR done previously
file_in_tsmr <- "tsmr_all_results_suppl.csv"
df_tsmr <- fread(file_in_tsmr, sep=",", data.table = F, nThread = 10,
                    colClasses = c(rep("character", 19))) %>%
            dplyr::select(Phecode, Method, Beta, SE, `P-value`) %>%
            filter(Method == "Inverse variance weighted")
df_tsmr$Beta <- as.numeric(df_tsmr$Beta)


bic.plot_list <- list()
effect.plot_list <- list()
for (idx in 1:nrow(df.cml)){
    print(idx)
    phecode <- df.cml[idx, ]$Phecode
    outcome <- df_mr_result[df_mr_result$Phecode == phecode, ]$Outcome

    cML_result <- readRDS(paste0(dir_res, "cML_result.exp_MetSnoUKB.out_", phecode, ".rds"))

    df1 <- data.frame(number_of_invalid_snp = 0:(length(cML_result$BIC_vec) - 1),
                    BIC = cML_result$BIC_vec)
    n_invalid_snp <- length(cML_result$BIC_invalid)
    df1$status <- factor(df1$number_of_invalid_snp == n_invalid_snp, levels = c(TRUE, FALSE))

    ### BIC scatter plot
    p1 <- ggplot(data = df1, aes(x = number_of_invalid_snp, y = BIC)) +
        geom_point(size = 1.5, color = "black") +
        geom_point(data = df1 %>% filter(status == TRUE), 
                   aes(x = number_of_invalid_snp, y = BIC),
                   color = "red")+
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        theme_light() +
        labs(x = "Number of invalid SNP",
            y = "BIC",
            title = paste0("Phecode: ", phecode, "; ", outcome),
            subtitle = paste0("N invalid IV = ", n_invalid_snp)) +
        theme(axis.text = element_text(size = 15),
                axis.title = element_text(size = 15, face = "bold"),
                plot.title = element_text(size = 18, face = "bold"),
                plot.subtitle = element_text(size = 14),
                legend.position = "None",
              panel.border = element_rect(color = "black"))
    
    
    ### Causal effect plot
    # Two-sample MR IVW
    df_plot <- readRDS(paste0(dir_cml_input, "exp_MetSnoUKB.out_", phecode, ".rds"))
    df_tsmr_phecode <- df_tsmr[df_tsmr$Phecode == phecode, ]
    slope_cml <- ifelse(df.cml[idx, ]$approach_selection == "MA_BIC", 
                        df.cml[idx, ]$MA_BIC_theta, df.cml[idx, ]$MA_BIC_DP_theta)
    

    if (df.cml[idx, ]$BIC_invalid_N == 0){ # No invalid SNP
        p2 <- ggplot(df_plot, aes(x = b_exp, y = b_out)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
            geom_errorbar(aes(xmin = b_exp-(1.96*se_exp), xmax = b_exp+(1.96*se_exp)),
                          color = "gray66", alpha = 0.5) +
            geom_errorbar(aes(ymin = b_out-(1.96*se_out), ymax = b_out+(1.96*se_out)),
                          color = "gray66", alpha = 0.5) +
            geom_point(color = "black", size = 1.5) +
            scale_x_continuous(limits = symmetric_limits) +
            scale_y_continuous(limits = symmetric_limits) +
            geom_abline(intercept = 0, slope = slope_cml, color = "red") +
            theme_light() +
            labs(x = "Beta exposure",
                y = "Beta outcome") +
            theme(axis.title = element_text(size = 15, face = "bold"),
                    axis.text = element_text(size = 15),
                    legend.position = "none",
                  panel.border = element_rect(color = "black")
                    )
    } else{
        # Plot dataframe
        validity <- c()
        for (idx2 in 1:nrow(df_plot)){
            valid <- TRUE
            if (idx2 %in% cML_result$BIC_invalid){
                valid <- FALSE
            }
            validity <- c(validity, valid)
        }
        df_plot$Validity <- factor(validity, levels = c(TRUE, FALSE))
        df_plot$tmp_group <- TRUE

        # Plot
        p2 <- ggplot(df_plot, aes(x = b_exp, y = b_out,
                    color = Validity)) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
                geom_errorbar(aes(xmin = b_exp-(1.96*se_exp), xmax = b_exp+(1.96*se_exp)),
                              color = "gray66", alpha = 0.5) +
                geom_errorbar(aes(ymin = b_out-(1.96*se_out), ymax = b_out+(1.96*se_out)),
                              color = "gray66", alpha = 0.5) +
                geom_point(size = 1.5) +
                scale_color_manual(values = c("black", "green")) +
                scale_x_continuous(limits = symmetric_limits) +
                scale_y_continuous(limits = symmetric_limits) +
            geom_abline(intercept = 0, slope = df_tsmr_phecode$Beta, color = "blue") +
                geom_abline(intercept = 0, slope = slope_cml, color = "red") +
                theme_light() +
                labs(x = "Beta exposure",
                    y = "Beta outcome") +
              theme(axis.title = element_text(size = 15, face = "bold"),
                    axis.text = element_text(size = 15),
                    legend.position = "none",
                    panel.border = element_rect(color = "black")
              )
    }
    
    ### Store two plot
    bic.plot_list[[idx]] <- p1
    effect.plot_list[[idx]] <- p2
}


#### Merge plots
combined1 <- bic.plot_list[[1]] + effect.plot_list[[1]] +
    bic.plot_list[[2]] + effect.plot_list[[2]] +
    plot_layout(nrow = 2, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot1.pdf"), 
       combined1, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 160
       )

combined2 <- bic.plot_list[[3]] + effect.plot_list[[3]] +
    bic.plot_list[[4]] + effect.plot_list[[4]] +
    bic.plot_list[[5]] + effect.plot_list[[5]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot2.pdf"), 
       combined2, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined3 <- bic.plot_list[[6]] + effect.plot_list[[6]] +
    bic.plot_list[[7]] + effect.plot_list[[7]] +
    bic.plot_list[[8]] + effect.plot_list[[8]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot3.pdf"), 
       combined3, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined4 <- bic.plot_list[[9]] + effect.plot_list[[9]] +
    bic.plot_list[[10]] + effect.plot_list[[10]] +
    bic.plot_list[[11]] + effect.plot_list[[11]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot4.pdf"), 
       combined4, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined5 <- bic.plot_list[[12]] + effect.plot_list[[12]] +
    bic.plot_list[[13]] + effect.plot_list[[13]] +
    bic.plot_list[[14]] + effect.plot_list[[14]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot5.pdf"), 
       combined5, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined6 <- bic.plot_list[[15]] + effect.plot_list[[15]] +
    bic.plot_list[[16]] + effect.plot_list[[16]] +
    bic.plot_list[[17]] + effect.plot_list[[17]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot6.pdf"), 
       combined6, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined7 <- bic.plot_list[[18]] + effect.plot_list[[18]] +
    bic.plot_list[[19]] + effect.plot_list[[19]] +
    bic.plot_list[[20]] + effect.plot_list[[20]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot7.pdf"), 
       combined7, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined8 <- bic.plot_list[[21]] + effect.plot_list[[21]] +
    bic.plot_list[[22]] + effect.plot_list[[22]] +
    bic.plot_list[[23]] + effect.plot_list[[23]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot8.pdf"), 
       combined8, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined9 <- bic.plot_list[[24]] + effect.plot_list[[24]] +
    bic.plot_list[[25]] + effect.plot_list[[25]] +
    bic.plot_list[[26]] + effect.plot_list[[26]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot9.pdf"), 
       combined9, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

combined10 <- bic.plot_list[[27]] + effect.plot_list[[27]] +
    bic.plot_list[[28]] + effect.plot_list[[28]] +
    bic.plot_list[[29]] + effect.plot_list[[29]] +
    plot_layout(nrow = 3, ncol = 2, byrow = TRUE)

ggsave(paste0(dir_out, "combined_plot10.pdf"), 
       combined10, device = "pdf", scale = 2,
       units = "mm", width = 180, height = 240
)

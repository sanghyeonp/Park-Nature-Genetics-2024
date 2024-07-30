
dir_MetS_PRS <- "../../PRS-CS_UKBB_target/"
dir_FOF_PRS <- "../../PRS-CS_UKBB_target/firstorder/"
basedir_IP_PRS <- "../../PRS-CS_UKBB_target/individual_phenotypes/"

PRS_files <- c(
    paste0(basedir_IP_PRS, 'BMI/', 'score_sum_indiv.profile.reformat'),
    paste0(basedir_IP_PRS, 'FG/', 'score_sum_indiv.profile.reformat'),
    paste0(basedir_IP_PRS, 'HDL/', 'score_sum_indiv.profile.reformat'),
    paste0(basedir_IP_PRS, 'HTN/', 'score_sum_indiv.profile.reformat'),
    paste0(basedir_IP_PRS, 'T2D/', 'score_sum_indiv.profile.reformat'),
    paste0(basedir_IP_PRS, 'TG/', 'score_sum_indiv.profile.reformat'),
    paste0(basedir_IP_PRS, 'WC/', 'score_sum_indiv.profile.reformat'),
    paste0(dir_FOF_PRS, 'score_sum_indiv_F1.profile.reformat'),
    paste0(dir_FOF_PRS, 'score_sum_indiv_F2.profile.reformat'),
    paste0(dir_FOF_PRS, 'score_sum_indiv_F3.profile.reformat'),
    paste0(dir_MetS_PRS, 'score_sum_indiv.profile.reformat')
)

PRS_names <- c("PRS_BMI", "PRS_FG", "PRS_HDL", "PRS_HTN", "PRS_T2D", "PRS_TG", "PRS_WC", 
                "PRS_F1", "PRS_F2", "PRS_F3", "PRS_MetS"
                )


read_calculated_PRS <- function(){
    cat("Reading PRS...\n")
    df <- NULL
    for (i in 1:length(PRS_names)){
        df_temp <- read.csv(PRS_files[i], sep='\t', header=TRUE)
        colnames(df_temp) <- c('FID', 'eid', 'PHENO', 'CNT', 'CNT2', 'PRS_raw')
        df_temp <- df_temp[, c("eid", 'PRS_raw')]
        # Normalize PRS
        df_temp$PRS_raw <- scale(df_temp$PRS_raw)
        colnames(df_temp) <- c('eid', PRS_names[i])

        if (i == 1){
            df <- df_temp
        } else{
            df <- merge(df, df_temp, by.x='eid', by.y='eid', all.x=TRUE, all.y=TRUE)
        }
    }
    return (df)
}

covariate_file <- "UKB_covariates.csv"

read_covariates <- function(){
    cat("Reading covariates...\n")
    df <- read.csv(covariate_file, header = TRUE, sep='\t')
    return (df)
}


merge_prs_covariates <- function(prs_df, covar_df){
    cat("Merging PRS and covariates...\n")
    df <- merge(prs_df, covar_df, by.x='eid', by.y='f.eid', all.x=TRUE, all.y=TRUE)
    return (df)
}

mets_clinical_file <- "UKB_MetS_clinical_phenotype.csv"

read_clinical_MetS <- function(){
    cat("Reading clinical MetS phenotype...\n")
    df <- read.csv(mets_clinical_file, header = TRUE)
    return (df)
}


merge_clinical_MetS <- function(df, cli_MetS_df){
    cat("Merging MetS clinical phenotype...\n")
    df <- merge(df, cli_MetS_df, by.x='eid', by.y='f.eid', all.x=TRUE, all.y=TRUE)
    return (df)
}



load_all_data <- function(){
    prs.df <- read_calculated_PRS()
    covar.df <- read_covariates()
    df <- merge_prs_covariates(prs.df, covar.df)
    clinic.mets.df <- read_clinical_MetS()
    df <- merge_clinical_MetS(df, clinic.mets.df)
    return (df)
}
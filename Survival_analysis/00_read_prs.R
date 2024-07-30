library(data.table)
library(dplyr)

read_PRS <- function(){
    ### MetS PRS
    file_in.prs_MetS <- "PRS_MetS_noUKB.csv"
    df_prs.MetS <- fread(file_in.prs_MetS, sep=",", data.table = F, nThread = 5) %>%
        mutate(FID = eid,
            IID = eid,
                PRS_scaled.MetS = scale(PRS),
            percentile.MetS = ntile(PRS_scaled.MetS, 100),
            PRS_group.MetS = ifelse(percentile.MetS %in% 1:20, 1, 
                            ifelse(percentile.MetS %in% 21:80, 2,
                                    ifelse(percentile.MetS %in% 81:99, 3, 
                                            ifelse(percentile.MetS == 100, 4, NA))))) %>%
        rename(PRS.MetS = PRS) %>%
        dplyr::select(FID, IID, PRS_scaled.MetS, percentile.MetS, PRS_group.MetS)


    ### BMI PRS
    file_in.prs_BMI <- "./BMI_score_sum_indiv.profile.reformat"
    df_prs.BMI <- fread(file_in.prs_BMI, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.BMI = SCORESUM) %>%
        mutate(PRS_scaled.BMI = scale(PRS.BMI),
            percentile.BMI = ntile(PRS_scaled.BMI, 100),
            PRS_group.BMI = ifelse(percentile.BMI %in% 1:20, 1, 
                                    ifelse(percentile.BMI %in% 21:80, 2,
                                            ifelse(percentile.BMI %in% 81:99, 3, 
                                                    ifelse(percentile.BMI == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.BMI, percentile.BMI, PRS_group.BMI)

    ### WC PRS
    file_in.prs_WC <- "./WC_score_sum_indiv.profile.reformat"
    df_prs.WC <- fread(file_in.prs_WC, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.WC = SCORESUM) %>%
        mutate(PRS_scaled.WC = scale(PRS.WC),
            percentile.WC = ntile(PRS_scaled.WC, 100),
            PRS_group.WC = ifelse(percentile.WC %in% 1:20, 1, 
                                    ifelse(percentile.WC %in% 21:80, 2,
                                            ifelse(percentile.WC %in% 81:99, 3, 
                                                    ifelse(percentile.WC == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.WC, percentile.WC, PRS_group.WC)

    ### HTN PRS
    file_in.prs_HTN <- "./HTN_score_sum_indiv.profile.reformat"
    df_prs.HTN <- fread(file_in.prs_HTN, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.HTN = SCORESUM) %>%
        mutate(PRS_scaled.HTN = scale(PRS.HTN),
            percentile.HTN = ntile(PRS_scaled.HTN, 100),
            PRS_group.HTN = ifelse(percentile.HTN %in% 1:20, 1, 
                                    ifelse(percentile.HTN %in% 21:80, 2,
                                            ifelse(percentile.HTN %in% 81:99, 3, 
                                                    ifelse(percentile.HTN == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.HTN, percentile.HTN, PRS_group.HTN)

    ### FG PRS
    file_in.prs_FG <- "./FG_score_sum_indiv.profile.reformat"
    df_prs.FG <- fread(file_in.prs_FG, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.FG = SCORESUM) %>%
        mutate(PRS_scaled.FG = scale(PRS.FG),
            percentile.FG = ntile(PRS_scaled.FG, 100),
            PRS_group.FG = ifelse(percentile.FG %in% 1:20, 1, 
                                    ifelse(percentile.FG %in% 21:80, 2,
                                            ifelse(percentile.FG %in% 81:99, 3, 
                                                    ifelse(percentile.FG == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.FG, percentile.FG, PRS_group.FG)

    ### T2D PRS
    file_in.prs_T2D <- "./T2D_score_sum_indiv.profile.reformat"
    df_prs.T2D <- fread(file_in.prs_T2D, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.T2D = SCORESUM) %>%
        mutate(PRS_scaled.T2D = scale(PRS.T2D),
            percentile.T2D = ntile(PRS_scaled.T2D, 100),
            PRS_group.T2D = ifelse(percentile.T2D %in% 1:20, 1, 
                                    ifelse(percentile.T2D %in% 21:80, 2,
                                            ifelse(percentile.T2D %in% 81:99, 3, 
                                                    ifelse(percentile.T2D == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.T2D, percentile.T2D, PRS_group.T2D)

    ### HDL PRS
    file_in.prs_HDL <- "./HDL_score_sum_indiv.profile.reformat"
    df_prs.HDL <- fread(file_in.prs_HDL, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.HDL = SCORESUM) %>%
        mutate(PRS_scaled.HDL = scale(PRS.HDL),
            percentile.HDL = ntile(PRS_scaled.HDL, 100),
            PRS_group.HDL = ifelse(percentile.HDL %in% 1:20, 1, 
                                    ifelse(percentile.HDL %in% 21:80, 2,
                                            ifelse(percentile.HDL %in% 81:99, 3, 
                                                    ifelse(percentile.HDL == 100, 4, NA)))),
            PRS_scaled.HDL.rev = -PRS_scaled.HDL,
            percentile.HDL.rev = ntile(PRS_scaled.HDL.rev, 100),
            PRS_group.HDL.rev = ifelse(percentile.HDL.rev %in% 1:20, 1, 
                                    ifelse(percentile.HDL.rev %in% 21:80, 2,
                                            ifelse(percentile.HDL.rev %in% 81:99, 3, 
                                                    ifelse(percentile.HDL.rev == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.HDL, percentile.HDL, PRS_group.HDL, 
                                PRS_scaled.HDL.rev, percentile.HDL.rev, PRS_group.HDL.rev)

    ### TG PRS
    file_in.prs_TG <- "./TG_score_sum_indiv.profile.reformat"
    df_prs.TG <- fread(file_in.prs_TG, sep="\t", data.table = F, nThread = 5) %>%
        dplyr::select(FID, IID, SCORESUM) %>%
        rename(PRS.TG = SCORESUM) %>%
        mutate(PRS_scaled.TG = scale(PRS.TG),
            percentile.TG = ntile(PRS_scaled.TG, 100),
            PRS_group.TG = ifelse(percentile.TG %in% 1:20, 1, 
                                    ifelse(percentile.TG %in% 21:80, 2,
                                            ifelse(percentile.TG %in% 81:99, 3, 
                                                    ifelse(percentile.TG == 100, 4, NA))))) %>%
        dplyr::select(FID, IID, PRS_scaled.TG, percentile.TG, PRS_group.TG)

    ### Merge
    df_prs <- merge(df_prs.MetS, df_prs.BMI, by = c("FID", "IID"), all.x = T)
    df_prs <- merge(df_prs, df_prs.WC, by = c("FID", "IID"), all.x = T)
    df_prs <- merge(df_prs, df_prs.HTN, by = c("FID", "IID"), all.x = T)
    df_prs <- merge(df_prs, df_prs.FG, by = c("FID", "IID"), all.x = T)
    df_prs <- merge(df_prs, df_prs.T2D, by = c("FID", "IID"), all.x = T)
    df_prs <- merge(df_prs, df_prs.HDL, by = c("FID", "IID"), all.x = T)
    df_prs <- merge(df_prs, df_prs.TG, by = c("FID", "IID"), all.x = T)

    return (df_prs)
}

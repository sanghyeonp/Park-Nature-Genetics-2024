library(data.table)
library(dplyr)

file_pheco_cov <- "KoGES_MetS_pheno_covar_table.txt"

df <- fread(file_pheco_cov, sep=" ", data.table=F, nThread=1) %>%
    mutate(MetS = ifelse(MetS == 1, 0, 1),
            MetS = ifelse(MetS == 2, 1, 0))

write.table(dplyr::select(df, FID, IID, MetS),
            file = "koges.pheno.txt",
            sep=" ", row.names=F, quote=F)

write.table(dplyr::select(df, -MetS),
            file = "koges.cov.txt",
            sep=" ", row.names=F, quote=F)

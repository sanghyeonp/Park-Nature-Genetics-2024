library(data.table)
library(dplyr)

### MetS
trait <- "MetS"
file_gwas <- "MetS_noUKB_GWAS.csv"

df <- fread(file_gwas, sep=",", data.table=F, nThread = 5) %>%
    dplyr::select(SNP, CHR, BP, A1, A2, MAF, est, Pval_Estimate) %>%
    rename(BETA = est,
           P = Pval_Estimate)

write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### TG
trait <- "TG"
file_gwas <- "CLEANED.TG_GLGC_noUKB.txt"

df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)

write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### HDL
trait <- "HDL"
file_gwas <- "CLEANED.HDL_GLGC_noUKB.txt"

df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)
        
write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### HTN
trait <- "HTN"
file_gwas <- "CLEANED.HTN_Finngen_r7.txt"

df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)
        
write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### FG
trait <- "FG"
file_gwas <- "CLEANED.FG_MAGIC.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)
        
write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### T2D
trait <- "T2D"
file_gwas <- "CLEANED.T2D_meta_FinngenR7_MVP.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)
        
write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### WC
trait <- "WC"
file_gwas <- "CLEANED.WC_GIANT2015.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)
        
write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)

### BMI
trait <- "BMI"
file_gwas <- "CLEANED.BMI_GIANT2015.txt"
df <- fread(file_gwas, sep="\t", data.table=F, nThread=10) %>%
    dplyr::select(SNP, CHR, POS, A1, A2, MAF, Effect, SE, Pval) %>%
    rename(BP = POS,
        BETA = Effect,
        P = Pval)
        
write.table(df, 
            file = paste0("gwas.", trait, "_noUKB.PRSice_input.txt"),
            sep=" ", row.names=F, quote=F)
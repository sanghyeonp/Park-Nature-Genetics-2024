require(GenomicSEM)
library(data.table)

NAME <- "MetS"

log.dir <- "./logs/"
dir.create(log.dir, showWarnings = FALSE)
log.dir <- paste0(log.dir, NAME)
dir.create(file.path(log.dir), showWarnings = FALSE)
setwd(log.dir)

sumstat.dir <- "../../gwas/EasyQC/results/"
reference.file <- "reference.1000G.maf.0.005.txt"

save.dir <- "./"

BMI_GIANT18 <- paste0(sumstat.dir, "CLEANED.BMI_GIANT2018.txt")
WC_Watanabe_UKB <- paste0(sumstat.dir, "CLEANED.WC_UKB.txt")
HDL_rev_GLGC <- paste0(sumstat.dir, "CLEANED.HDL_GLGC_UKB.txt")
TG_GLGC <- paste0(sumstat.dir, "CLEANED.TG_GLGC_UKB.txt")
FG_Lagou_MAGIC_noUKB <- paste0(sumstat.dir, "CLEANED.FG_MAGIC.txt")
T2D_meta_FinR7_Mahajan22_MVP <- paste0(sumstat.dir, "CLEANED.T2D_meta_FinngenR7_Mahajan2022_MVP.txt")
HYP_meta_FinR7_UKB <- paste0(sumstat.dir, "CLEANED.HTN_meta_FinngenR7_UKB.txt")

files=c(BMI_GIANT18, WC_Watanabe_UKB, HDL_rev_GLGC, TG_GLGC, FG_Lagou_MAGIC_noUKB,
        T2D_meta_FinR7_Mahajan22_MVP, HYP_meta_FinR7_UKB)
trait.names=c("BMI_GIANT18", "WC_Watanabe_UKB", "HDL_rev_GLGC", "TG_GLGC", "FG_Lagou_MAGIC_noUKB",
        "T2D_meta_FinR7_Mahajan22_MVP", "HYP_meta_FinR7_UKB")

se.logit=c(F, F, F, F, F, F, F)
OLS=c(T, T, T, T, T, T, T)
linprob=c(F, F, F, F, F, F, F)

N=c(NA, 385932, NA, NA, NA, 597437, 508612)
info.filter=0.6
maf.filter=0.01
keep.indel=FALSE
parallel=TRUE
cores=7

MetS_sumstat <-sumstats(files=files,
                        ref=reference.file, 
                        trait.names=trait.names, 
                        se.logit=se.logit, OLS=OLS, linprob=linprob, N=N, 
                        info.filter=info.filter, maf.filter=maf.filter, 
                        keep.indel=keep.indel, parallel=parallel, cores=cores
                        )

save(MetS_sumstat, file=paste(save.dir, NAME, "_sumstat.RData", sep=""))

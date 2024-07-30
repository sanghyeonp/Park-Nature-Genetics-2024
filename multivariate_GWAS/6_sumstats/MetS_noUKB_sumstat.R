require(GenomicSEM)
library(data.table)

NAME <- "MetS_noUKB"

log.dir <- "./logs/"
dir.create(log.dir, showWarnings = FALSE)
log.dir <- paste0(log.dir, NAME)
dir.create(file.path(log.dir), showWarnings = FALSE)
setwd(log.dir) # Place where the logs will be saved

sumstat.dir <- "../../gwas/EasyQC/results/"
reference.file <- "reference.1000G.maf.0.005.txt"

save.dir <- "./"

BMI_GIANT15_noUKB <- paste0(sumstat.dir, "CLEANED.BMI_GIANT2015.txt")
WC_GIANT_noUKB <- paste0(sumstat.dir, "CLEANED.WC_GIANT2015.txt")
HDL_rev_GLGC_noUKB <- paste0(sumstat.dir, "CLEANED.HDL_GLGC_noUKB.txt")
TG_GLGC_noUKB <- paste0(sumstat.dir, "CLEANED.TG_GLGC_noUKB.txt")
FG_Lagou_MAGIC_noUKB <- paste0(sumstat.dir, "CLEANED.FG_MAGIC.txt")
T2D_meta_FinR7_MVP_noUKB <- paste0(sumstat.dir, "CLEANED.T2D_meta_FinngenR7_MVP.txt")
HYP_FinR7_noUKB <- paste0(sumstat.dir, "CLEANED.T2D_Finngen_r7.txt")

files=c(BMI_GIANT15_noUKB, WC_GIANT_noUKB, HDL_rev_GLGC_noUKB, TG_GLGC_noUKB, FG_Lagou_MAGIC_noUKB,
        T2D_meta_FinR7_MVP_noUKB, HYP_FinR7_noUKB)
trait.names=c("BMI_GIANT15_noUKB", "WC_GIANT_noUKB", "HDL_rev_GLGC_noUKB", "TG_GLGC_noUKB", "FG_Lagou_MAGIC_noUKB",
        "T2D_meta_FinR7_MVP_noUKB", "HYP_FinR7_noUKB")

se.logit=c(F, F, F, F, F, F, T)
OLS=c(T, T, T, T, T, T, F)
linprob=c(F, F, F, F, F, F, F)

N=c(NA, NA, NA, NA, NA, 345698, 247289)
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

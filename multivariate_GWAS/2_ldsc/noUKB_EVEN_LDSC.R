require(GenomicSEM)
library(data.table)

NAME <- "MetS_noUKB_EVEN"

munge.dir <- "../1_munge/"
ldsc.output.dir <- "./ldscout/"
ld_weights <- "./eur_w_ld_chr/"
logs.dir <- "./logs/"
dir.create(file.path(logs.dir), showWarnings = FALSE)


BMI_GIANT15_noUKB <- paste0(munge.dir, "BMI_GIANT15_noUKB.sumstats.gz")
WC_GIANT_noUKB <- paste0(munge.dir, "WC_GIANT_noUKB.sumstats.gz")
HDL_rev_GLGC_noUKB <- paste0(munge.dir, "HDL_rev_GLGC_noUKB.sumstats.gz")
TG_GLGC_noUKB <- paste0(munge.dir, "TG_GLGC_noUKB.sumstats.gz")
FG_Lagou_MAGIC_noUKB <- paste0(munge.dir, "FG_Lagou_MAGIC_noUKB.sumstats.gz")
T2D_meta_FinR7_MVP_noUKB <- paste0(munge.dir, "T2D_meta_FinR7_MVP_noUKB.sumstats.gz")
HYP_FinR7_noUKB <- paste0(munge.dir, "HYP_FinR7_noUKB.sumstats.gz")

TS <- c(BMI_GIANT15_noUKB, WC_GIANT_noUKB, HDL_rev_GLGC_noUKB, TG_GLGC_noUKB, FG_Lagou_MAGIC_noUKB,
        T2D_meta_FinR7_MVP_noUKB, HYP_FinR7_noUKB)
sample.prev <- c(NA, NA, NA, NA, NA, 0.5, 0.5)
population.prev <- c(NA, NA, NA, NA, NA, 0.07, 0.442)

trait.names <- c("BMI_GIANT15_noUKB", "WC_GIANT_noUKB", "HDL_rev_GLGC_noUKB", "TG_GLGC_noUKB", "FG_Lagou_MAGIC_noUKB",
        "T2D_meta_FinR7_MVP_noUKB", "HYP_FinR7_noUKB")

MetS_ldscout <- ldsc(TS, sample.prev, population.prev, ld_weights, ld_weights, trait.names,
                        ldsc.log=paste0(logs.dir, NAME), stand=TRUE, select='EVEN')


save(MetS_ldscout, file=paste0(ldsc.output.dir, NAME, "_ldscout.RData"))

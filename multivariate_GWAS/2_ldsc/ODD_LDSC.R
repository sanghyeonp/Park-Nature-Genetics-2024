require(GenomicSEM)
library(data.table)

NAME <- "MetS_ODD"

munge.dir <- "../1_munge/"
ldsc.output.dir <- "./ldscout/"
ld_weights <- "./eur_w_ld_chr/"
logs.dir <- "./logs/"
dir.create(file.path(logs.dir), showWarnings = FALSE)


BMI_GIANT18 <- paste0(munge.dir, "BMI_GIANT18.sumstats.gz")
WC_Watanabe_UKB <- paste0(munge.dir, "WC_Watanabe_UKB.sumstats.gz")
HDL_rev_GLGC <- paste0(munge.dir, "HDL_rev_GLGC.sumstats.gz")
TG_GLGC <- paste0(munge.dir, "TG_GLGC.sumstats.gz")
FG_Lagou_MAGIC_noUKB <- paste0(munge.dir, "FG_Lagou_MAGIC_noUKB.sumstats.gz")
T2D_meta_FinR7_Mahajan22_MVP <- paste0(munge.dir, "T2D_meta_FinR7_Mahajan22_MVP.sumstats.gz")
HYP_meta_FinR7_UKB <- paste0(munge.dir, "HYP_meta_FinR7_UKB.sumstats.gz")



TS <- c(BMI_GIANT18, WC_Watanabe_UKB, HDL_rev_GLGC, TG_GLGC, FG_Lagou_MAGIC_noUKB,
        T2D_meta_FinR7_Mahajan22_MVP, HYP_meta_FinR7_UKB)   # Array of phenotypes
sample.prev <- c(NA, NA, NA, NA, NA, 0.5, 0.5)
population.prev <- c(NA, NA, NA, NA, NA, 0.07, 0.442)

trait.names <- c("BMI_GIANT18", "WC_Watanabe_UKB", "HDL_rev_GLGC", "TG_GLGC", "FG_Lagou_MAGIC_noUKB",
        "T2D_meta_FinR7_Mahajan22_MVP", "HYP_meta_FinR7_UKB")


MetS_ldscout <- ldsc(TS, sample.prev, population.prev, ld_weights, ld_weights, trait.names,
                        ldsc.log=paste0(logs.dir, NAME), stand=TRUE, select='ODD')

save(MetS_ldscout, file=paste0(ldsc.output.dir, NAME, "_ldscout.RData"))


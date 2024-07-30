require(GenomicSEM)
library(data.table)

sumstat.dir <- "../../gwas/EasyQC/results/"
munge.dir <- "./"
log.dir <- "./logs/"
dir.create(file.path(log.dir), showWarnings = FALSE)

reference_file <- 'w_hm3.snplist'

TG_GLGC <- paste0(sumstat.dir, "CLEANED.TG_GLGC_UKB.txt")
TG_GLGC_noUKB <- paste0(sumstat.dir, "CLEANED.TG_GLGC_noUKB.txt")

files.TG <- c(TG_GLGC, TG_GLGC_noUKB)
trait.names.TG <- c("TG_GLGC", "TG_GLGC_noUKB")
N.TG <- c(NA, NA)

files <- files.TG
trait.names <- trait.names.TG
N <- N.TG

HDL_rev_GLGC <- paste0(sumstat.dir, "CLEANED.HDL_GLGC_UKB.txt")
HDL_rev_GLGC_noUKB <- paste0(sumstat.dir, "CLEANED.HDL_GLGC_noUKB.txt")

files.HDL <- c(HDL_rev_GLGC, HDL_rev_GLGC_noUKB)
trait.names.HDL <- c("HDL_rev_GLGC", "HDL_rev_GLGC_noUKB")
N.HDL <- c(NA, NA)

files <- c(files, files.HDL)
trait.names <- c(trait.names, trait.names.HDL)
N <- c(N, N.HDL)

HYP_meta_FinR7_UKB <- paste0(sumstat.dir, "CLEANED.HTN_meta_FinngenR7_UKB.txt")
HYP_FinR7_noUKB <- paste0(sumstat.dir, "CLEANED.T2D_Finngen_r7_GC.txt")

files.HYP <- c(HYP_meta_FinR7_UKB, HYP_FinR7_noUKB)
trait.names.HYP <- c("HYP_meta_FinR7_UKB", "HYP_FinR7_noUKB")
N.HYP <- c(508612, 247289)

files <- c(files, files.HYP)
trait.names <- c(trait.names, trait.names.HYP)
N <- c(N, N.HYP)

FG_Lagou_MAGIC_noUKB <- paste0(sumstat.dir, "CLEANED.FG_MAGIC_A1asEffect.txt")

files.FG <- c(FG_Lagou_MAGIC_noUKB)
trait.names.FG <- c("FG_Lagou_MAGIC_noUKB")
N.FG <- c(NA)

files <- c(files, files.FG)
trait.names <- c(trait.names, trait.names.FG)
N <- c(N, N.FG)

T2D_meta_FinR7_Mahajan22_MVP <- paste0(sumstat.dir, "CLEANED.T2D_meta_FinngenR7_Mahajan2022_MVP.txt")
T2D_meta_FinR7_MVP_noUKB <- paste0(sumstat.dir, "CLEANED.T2D_meta_FinngenR7_MVP.txt")

files.T2D <- c(T2D_meta_FinR7_Mahajan22_MVP, T2D_meta_FinR7_MVP_noUKB)
trait.names.T2D <- c("T2D_meta_FinR7_Mahajan22_MVP", "T2D_meta_FinR7_MVP_noUKB")
N.T2D <- c(597437, 345698)

files <- c(files, files.T2D)
trait.names <- c(trait.names, trait.names.T2D)
N <- c(N, N.T2D)


WC_Watanabe_UKB <- paste0(sumstat.dir, "CLEANED.WC_UKB.txt")
WC_GIANT_noUKB <- paste0(sumstat.dir, "CLEANED.WC_GIANT2015.txt")

files.WC <- c(WC_Watanabe_UKB, WC_GIANT_noUKB)
trait.names.WC <- c("WC_Watanabe_UKB", "WC_GIANT_noUKB")
N.WC <- c(385932, NA)

files <- c(files, files.WC)
trait.names <- c(trait.names, trait.names.WC)
N <- c(N, N.WC)

BMI_GIANT18 <- paste0(sumstat.dir, "CLEANED.BMI_GIANT2018.txt")
BMI_GIANT15_noUKB <- paste0(sumstat.dir, "CLEANED.BMI_GIANT2015.txt")

files.BMI <- c(BMI_GIANT18, BMI_GIANT15_noUKB)
trait.names.BMI <- c("BMI_GIANT18", "BMI_GIANT15_noUKB")
N.BMI <- c(NA, NA)

files <- c(files, files.BMI)
trait.names <- c(trait.names, trait.names.BMI)
N <- c(N, N.BMI)

for (i in 1:length(files)){
    cat(paste0(files[i], "\n"))
    munge(  files=c(files[i]), 
            hm3 = reference_file, 
            trait.names=c(trait.names[i]), 
            N=c(N[i]), 
            info.filter = 0.9, 
            maf.filter = 0.01,
            parallel = FALSE,
            log.name=paste0(log.dir, trait.names[i])
        )
}

require(GenomicSEM)
library(data.table)

NAME <- "MetS_noUKB_EVEN"

model <- "F1 =~ NA*BMI_GIANT15_noUKB + WC_GIANT_noUKB
        F2 =~ NA*T2D_meta_FinR7_MVP_noUKB + FG_Lagou_MAGIC_noUKB + HYP_FinR7_noUKB
        F3 =~ NA*HDL_rev_GLGC_noUKB + TG_GLGC_noUKB
        MetS =~ NA*F1 + F2 + F3
        MetS ~~ 1*MetS
        F1 ~~ 1*F1
        F2 ~~ 1*F2
        F3 ~~ 1*F3
        F1 ~~ 0*F2
        F1 ~~ 0*F3
        F2 ~~ 0*F3
        "

ldsc.out <- "../2_ldsc/ldscout/MetS_noUKB_EVEN_ldscout.RData"

loadRData <- function(file){
    #loads an RData file, and returns it
    load(file)
    get(ls()[ls() != "file"])
}

ldscout <- loadRData(ldsc.out)
CFA <- usermodel(ldscout, estimation = "DWLS", model = model, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE)

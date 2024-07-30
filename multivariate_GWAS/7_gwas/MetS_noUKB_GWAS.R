require(GenomicSEM)
library(data.table)

NAME <- "MetS_noUKB"

ldsc.out <- "../2_ldsc/ldscout/MetS_noUKB_ldscout.RData"
sumstat.out <- "../7_sumstat/MetS_noUKB_sumstat.RData"

load(ldsc.out)
load(sumstat.out)

model <- "F1 =~ 1*BMI_GIANT15_noUKB + WC_GIANT_noUKB
        F2 =~ 1*T2D_meta_FinR7_MVP_noUKB + FG_Lagou_MAGIC_noUKB + HYP_FinR7_noUKB
        F3 =~ 1*HDL_rev_GLGC_noUKB + TG_GLGC_noUKB
        MetS =~ 1*F1 + F2 + F3
        F1 ~~ 0*F2
        F1 ~~ 0*F3
        F2 ~~ 0*F3
	T2D_meta_FinR7_MVP_noUKB ~~ a*T2D_meta_FinR7_MVP_noUKB
	WC_GIANT_noUKB ~~ b*WC_GIANT_noUKB
        BMI_GIANT15_noUKB ~~ c*BMI_GIANT15_noUKB
        HDL_rev_GLGC_noUKB ~~ d*HDL_rev_GLGC_noUKB
        HYP_FinR7_noUKB ~~ e*HYP_FinR7_noUKB
        FG_Lagou_MAGIC_noUKB ~~ f*FG_Lagou_MAGIC_noUKB
        TG_GLGC_noUKB ~~ g*TG_GLGC_noUKB
        F1 ~~ h*F1
        F2 ~~ i*F2
        F3 ~~ j*F3
        a > 0.001
        b > 0.001
        c > 0.001
        d > 0.001
        e > 0.001
        f > 0.001
        g > 0.001
        h > 0.001
        i > 0.001
        j > 0.001
        MetS ~ SNP
        "

MetS.GWAS <- userGWAS(covstruc=MetS_ldscout, 
                        SNPs=MetS_sumstat,
                        estimation="DWLS", 
                        model=model,
                        cores=16, 
                        toler=1e-50, 
                        SNPSE=FALSE, 
                        sub=c('MetS ~ SNP'),
                        parallel=TRUE, 
                        GC="standard", 
                        MPI=FALSE, 
                        smooth_check=FALSE,
                        printwarn=TRUE)

write.table(MetS.GWAS, 
        paste0(NAME, "_GWAS.tsv"),
        sep="\t", row.names=F, quote=F)

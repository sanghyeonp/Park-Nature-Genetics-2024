library(data.table)
library(dplyr)

trait_list <- c("MetS", "TG", "HDL",
                "FG", 
                "BMI",
                "HTN",
                "T2D", "WC",
                "F1", "F2", "F3",
                "MetS_lind", "MetS_walree")


for (trait in trait_list){
    print(trait)
    dir_cojo <- paste0("../cojo.", trait, "/")
    file_in <- paste0(dir_cojo, "cojo.", trait, ".jma.cojo")
    df <- fread(file_in, sep="\t", data.table= F, nThread = 1)

    write.table(df,
                paste0("cojo_out_jma.", trait, ".csv"),
                sep=",", row.names = F, quote = F
    )

    write.table(df$SNP,
                paste0("cojo_out_jma.", trait, ".snplist"),
                row.names = F, quote = F, col.names = F)

    write.table(df[, c("SNP", "Chr", "bp", "p")],
                paste0("cojo_out_jma.", trait, ".snplist.chr_pos"),
                row.names = F, quote = F)
}

library(data.table)
library(dplyr)


### MetS
trait <- "MetS"
file_in <- "leadSNPs.txt"
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### TG
trait <- "TG"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### HDL
trait <- "HDL"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### HTN
trait <- "HTN"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### FG
trait <- "FG"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### T2D
trait <- "T2D"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### WC
trait <- "WC"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)


### BMI
trait <- "BMI"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### MetS Qsnp
trait <- "MetS_Qsnp"
# file_in <- "./clumping.MetS_Qsnp/MetS_Qsnp.clumped"
file_in <- "./FUMA.MetS_Qsnp/leadSNPs_MetS12_Qsnp.txt"
df <- fread(file_in, data.table= F, nThread = 1)
write.table(df$SNP,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### F1
trait <- "F1"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### F2
trait <- "F2"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### F3
trait <- "F3"
file_in <- paste0("leadSNPs_", trait, ".txt")
df <- fread(file_in, data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)
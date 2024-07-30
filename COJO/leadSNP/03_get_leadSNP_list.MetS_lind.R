library(data.table)
library(dplyr)


### MetS
trait <- "MetS_lind"
file_in <- "leadSNP.MetS_lind.txt"
df <- fread(file_in, data.table= F, nThread = 1)
write.table(df$SNP,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)

### MetS
trait <- "MetS_walree"
file_in <- "leadSNP.MetS_walree.txt"
df <- fread(file_in, data.table= F, nThread = 1)
write.table(df$SNP,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)


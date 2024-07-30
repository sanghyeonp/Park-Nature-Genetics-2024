library(data.table)
library(dplyr)


### MetS
trait <- "MetS.KoGES"
file_in <- "leadSNPs.MetS.KoGES.txt"
df <- fread(file_in, sep="\t", data.table= F, nThread = 1)
write.table(df$rsID,
            paste0("leadSNP_list.", trait, ".txt"),
            sep="\t", row.names = F, quote = F, col.names = F
)


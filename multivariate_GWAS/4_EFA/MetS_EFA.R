require(GenomicSEM)
library(data.table)
require(Matrix)
require(stats)


NAME <- "MetS"
N_FACTOR <- 3


ldsc.out <- paste0("../2_ldsc/ldscout/", NAME, "_ldscout.RData")


loadRData <- function(file){
    load(file)
    get(ls()[ls() != "file"])
}

ldscout <- loadRData(ldsc.out)

ROTATION <- "promax"
EFA.log <- paste0(NAME, "_", ROTATION, "_EFA.log")
sink(EFA.log, split=TRUE) 
S <- ldscout$S

if(length(which(S<0, arr.ind = TRUE)) != 0){
    S <- as.matrix((nearPD(S, corr = FALSE))$mat)
    cat("S matrix has been smoothed!!!\n")
}


for(i in 1:N_FACTOR){
    EFA <- factanal(covmat=S, factors=i, rotation=ROTATION)
    print(EFA)
    write.csv(EFA$uniqueness, paste0(NAME, "_", ROTATION, "_EFA", i, "_uniqueness.csv"), row.names=T)
    write.csv(EFA$loadings, paste0(NAME, "_", ROTATION, "_EFA", i, "_loadings.csv"), row.names=T)
}

cat(paste0("Check log file: ", EFA.log, "\n"))
sink()

ROTATION <- "varimax"
EFA.log <- paste0(NAME, "_", ROTATION, "_EFA.log")
sink(EFA.log, split=TRUE)
S <- ldscout$S

if(length(which(S<0, arr.ind = TRUE)) != 0){
    S <- as.matrix((nearPD(S, corr = FALSE))$mat)
    cat("S matrix has been smoothed!!!\n")
}


for(i in 1:N_FACTOR){
    EFA <- factanal(covmat=S, factors=i, rotation=ROTATION)
    print(EFA)
    write.csv(EFA$uniqueness, paste0(NAME, "_", ROTATION, "_EFA", i, "_uniqueness.csv"), row.names=T)
    write.csv(EFA$loadings, paste0(NAME, "_", ROTATION, "_EFA", i, "_loadings.csv"), row.names=T)
}

cat(paste0("Check log file: ", EFA.log, "\n"))
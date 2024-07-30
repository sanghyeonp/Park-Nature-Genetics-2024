suppressPackageStartupMessages({
    require(argparse)
    library(EFAtools)
    library(data.table)
    library(stringr)
})

parser <- argparse::ArgumentParser(description=":: Computes Kaiser-Meyer-Olkin criterion of a Genomic SEM phenotype set ::", formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--name", required=FALSE, default="NA",
                    help="Name of the model.")
parser$add_argument("--rdata", required=TRUE, 
                    help="Path Rdata generated using ldsc() function in Genomic SEM.")
parser$add_argument("--outDir", required=TRUE, 
                    help="Path KMO result will be saved.")

## Get parser arguments
args <- parser$parse_args()

if(str_sub(args$outDir,-1,-1) != '/'){
    args$outDir <- paste0(args$outDir, '/')
}
if(args$name == 'NA'){
    args$name <- strsplit(basename(args$rdata), '[.]')[[1]][1]
}

loadRData <- function(file){
    #loads an RData file, and returns it
    load(file)
    get(ls()[ls() != "file"])
}

ldscout <- loadRData(file=args$rdata)

cor.matrix <- ldscout$S_Stand

kmo.result <- KMO(x=data.matrix(cor.matrix))

print(kmo.result)

df <- data.frame(kmo.result$KMO_i)
df.final <- transpose(df)
colnames(df.final) <- rownames(df)
rownames(df.final) <- colnames(df)
df.final <- cbind(Overall_KMO = kmo.result$KMO, df.final)
rownames(df.final) <- "KMO value"

write.csv(df.final, file=paste0(args$outDir, args$name, '.csv'))

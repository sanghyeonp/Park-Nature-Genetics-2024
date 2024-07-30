library(TwoSampleMR)
library(MRInstruments)

### Exposure data preparation
exp_clumped_file_path <- "MetS_exposure_clumped.rds"

if (file.exists(exp_clumped_file_path) == FALSE){
  exposure_file <- "MetS_exposure.csv"
  
  exp_dat <- read_exposure_data(
    filename = exposure_file,
    sep = ",",
    snp_col = "SNP",
    beta_col = "est",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "MAF",
    pval_col = "Pval_Estimate",
    units_col = "Units",
    samplesize_col = "N",
    chr_col = "CHR",
    pos_col = "BP"
  )
  
  
  exp_dat_clump <- clump_data(exp_dat,
                              clump_kb = 10000,
                              clump_r2 = 0.001,
                              clump_p1 = 5e-8,
                              clump_p2 = 5e-8,
                              pop = "EUR")
  
  
  saveRDS(exp_dat_clump, file = exp_clumped_file_path)
}


### Harmonized data preparation
exp_dat <- readRDS(exp_clumped_file_path)


outcome_savedir <- "./outcome_rds/"
harmonized_savedir <- "./harmonized_rds/"

meta_file <- "TSMR_metadata_for_analysis.tsv"

metadata <- read.table(meta_file, sep="\t", row.names=NULL, header=TRUE, quote="")

outcome_dir <- "./reformat_outcome_gwas/"

N.IV <- NULL
for (idx in 1:nrow(metadata)){
  print(idx)
  print(metadata$Phenotype.from.PRS.PheWAS[idx])
  print(metadata$File[idx])
  print("#############################")
  
  outcome_rdsfile <- paste0(outcome_savedir, metadata$File[idx], ".outcome.rds")
  harmonized_rdsfile <- paste0(harmonized_savedir, metadata$File[idx], ".harmonized.rds")
  
  if (file.exists(harmonized_rdsfile) == TRUE){
    dat0 <- readRDS(harmonized_rdsfile)
  } else{
    outcome_file <- paste0(outcome_dir, metadata$File[idx], ".tsv")
    
    outcome_dat <- read_outcome_data(
      snps = exp_dat$SNP,
      filename = outcome_file,
      sep = "\t",
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "minor_allele",
      other_allele_col = "major_allele",
      eaf_col = "minor_AF",
      pval_col = "pval",
      units_col = "Units"
    )
    
    saveRDS(outcome_dat, file=outcome_rdsfile)
    
    dat0 <- harmonise_data(
      exposure_dat = exp_dat, 
      outcome_dat = outcome_dat
    )
    
    saveRDS(dat0, file=harmonized_rdsfile)
  }
  N.IV <- c(N.IV, nrow(dat0[dat0$mr_keep == TRUE, ]))
}

metadata$N.IV <- N.IV

write.table(metadata, file="metadata_harmonized.tsv",
            sep="\t", row.names=FALSE)

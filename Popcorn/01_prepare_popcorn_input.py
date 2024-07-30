import pandas as pd

file = "MetS_GWAS_MetS.csv"
df_ukb = pd.read_csv(file, sep=",", index_col=False)
df_ukb = df_ukb[['SNP', 'CHR', 'BP', 'A1', 'A2', 'MAF', 'est', 'SE', 'Pval_Estimate']]
n_init = len(df_ukb)
df_ukb = df_ukb[df_ukb['MAF'] > 0.01]
print("Number of SNPs dropped by MAF filter (> 0.01): {:,}".format(n_init - len(df_ukb)))
df_ukb.drop(columns=['MAF'], inplace=True)
df_ukb.columns = ['SNP', 'chr', 'pos', 'A1', 'A2', 'beta', 'SE', 'P']
df_ukb['N'] = 1384348

df_ukb.to_csv("MetS_GWAS_MetS_maf0.01_popcornIN.tsv", sep="\t", index=False)

file = "KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid.modified"
df_koges = pd.read_csv(file, sep="\t", index_col=False)
df_koges = df_koges[['RSID', '#CHROM', 'POS', 'effect_allele', 'other_allele', 'MAF', 'BETA', 'LOG(OR)_SE', 'P']]
n_init = len(df_koges)
df_koges = df_koges[df_koges['MAF'] > 0.01]
print("Number of SNPs dropped by MAF filter (> 0.01): {:,}".format(n_init - len(df_koges)))
df_koges.drop(columns=['MAF'], inplace=True)
df_koges.columns = ['SNP', 'chr', 'pos', 'A1', 'A2', 'beta', 'SE', 'P']
df_koges['N'] = 62314

df_koges.to_csv("KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid.modified_maf0.01_popcornIN.tsv", sep="\t", index=False)

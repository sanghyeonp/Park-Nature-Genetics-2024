import pandas as pd
import code

df1 = pd.read_csv("KoGES2022_MetS_binarization_phenotype_table.csv", sep=",", index_col=False)
df1 = df1[['DIST_ID', 'MetS', 'AGE', 'SEX']]
df1.rename(columns={"DIST_ID":"IID"}, inplace=True)
df1['MetS'] = df1['MetS'].apply(lambda x: 'case' if x == 1 else 'control')
df1['MetS'] = df1['MetS'].apply(lambda x: 2 if x == 'case' else 1)

# PC file
df2 = pd.read_csv("KCHIP.abc.snp.pca.eigenvec", sep="\t", index_col=False)
df2.rename(columns={"#FID":"FID"}, inplace=True)

# Merge
df = df1.merge(df2, how="left", on="IID")

# Check missing value
df_na = df[df.isnull().any(axis=1)] # 86명 PC값이 없음. Sample QC에서 탈락된 사람들.

# Drop na
df.dropna(axis=0, how="any", inplace=True)

# Save
df = df[['FID', 'IID', 'MetS', 'AGE', 'SEX', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
    'PC6', 'PC7', 'PC8', 'PC9', 'PC10']]
df.to_csv("KoGES_MetS_pheno_covar_table.txt", sep=" ", index=False)


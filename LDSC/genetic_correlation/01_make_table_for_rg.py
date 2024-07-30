import pandas as pd

metadata_final_file = "metadata_final_forLDSC_20221229.csv"

df = pd.read_csv(metadata_final_file, sep=",", index_col=False)

mets_munge_file = "../MetS_GWAS_munge.sumstats.gz"

df2 = df[['Category', 'Trait']]
df2['MetSMunge'] = [mets_munge_file] * len(df)

df2['TraitMunge'] = df['Munged_path'].tolist()

df2['Trait'] = df2['Trait'].apply(lambda trait: trait.replace(",", "_"))

df2.to_csv("./metadata_for_LDSC_rg.csv", sep=",", index=False)

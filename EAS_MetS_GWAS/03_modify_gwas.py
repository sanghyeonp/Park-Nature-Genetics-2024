import pandas as pd
import math

annovar_annot_file = "KCHIP.snp.imputed.HRC.abc.anno"
df_annot = pd.read_csv(annovar_annot_file, sep="\t", index_col=False)
df_annot = df_annot[['SNP', 'avsnp150']]
df_annot.rename(columns={"avsnp150":"RSID"}, inplace=True)

plink_result = "KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid"
df = pd.read_csv(plink_result, sep="\t", index_col=False)

df = df.merge(df_annot, how="left", left_on="ID", right_on="SNP")

df['BETA'] = df['OR'].apply(lambda x: math.log(x))

maf_file="KCHIP.snp.imputed.HRC.abc.frq"

def read_formatted_file(file):
    with open(file, 'r') as f:
        rows = [row.strip() for row in f.readlines()]
        data = []
        for row in rows:
            row = row.split(sep=" ")
            row = [ele for ele in row if ele]
            data.append(row)
    
    df = pd.DataFrame(data[1:], columns=data[0])
    return df

df_maf = read_formatted_file(maf_file)
df_maf = df_maf[['SNP', 'MAF']]
df = df.merge(df_maf, how="left", left_on="ID", right_on="SNP")


def get_a1_a2(ref, alt, a1):
    if alt == a1:
        return [a1, ref]
    return [ref, a1]

df[['effect_allele', 'other_allele']] = df.apply(lambda row: get_a1_a2(row['REF'], row['ALT'], row['A1']), axis=1, result_type='expand')

## Save
df.to_csv("KoGES_MetS_logistic_plink.MetS.glm.logistic.hybrid.modified", sep="\t", index=False)

import os
import pandas as pd
import subprocess
from multiprocessing import Pool, Manager


metadata_file = "../01_metadata_prep/TSMR_metadata_for_analysis.tsv"
df_metadata = pd.read_csv(metadata_file, sep="\t", index_col=False)

ukb_gwas_dir = "./GWASsumstat/UKB"

savedir = "./reformat_outcome_gwas/"

def read_variant_annotation_file():
    print("Reading variant annotation file...")
    df_annotation = pd.read_csv("variants.tsv", sep="\t", index_col=False)  
    df_annotation = df_annotation[['variant', 'rsid']]
    return df_annotation


def reformat_ukb_for_mr(sumstat, annotation_df, savedir):
    _, file = os.path.split(sumstat)
    
    if os.path.exists(os.path.join(savedir, file)):
        return None

    print("Reformatting {}".format(sumstat))
    df = pd.read_csv(sumstat, sep="\t", index_col=False)
    df = df[['variant', 'minor_allele', 'minor_AF', 'n_complete_samples', 'beta', 'se', 'pval']]

    df['CHR'] = df['variant'].apply(lambda x: x.split(sep=":")[0])
    df['POS'] = df['variant'].apply(lambda x: x.split(sep=":")[1])
    df['A1'] = df['variant'].apply(lambda x: x.split(sep=":")[2])
    df['A2'] = df['variant'].apply(lambda x: x.split(sep=":")[3])

    df['A1 is MA'] = df.apply(lambda row: row['A1'] == row['minor_allele'], axis=1)
    df['A1 new'] = df.apply(lambda row: df['A1'] if row['A1 is MA'] else row['A2'], axis=1)
    df['A2 new'] = df.apply(lambda row: df['A2'] if row['A1 is MA'] else row['A1'], axis=1)


    df2 = pd.merge(df, annotation_df, on='variant', how="left")

    df2 = df2[['rsid', 'CHR', 'POS', 'A1 new', 'A2 new', 'minor_AF', 'beta', 'se', 'pval', 'n_complete_samples']]
    df2.columns = ['rsid', 'chr', 'pos', 'minor_allele', 'major_allele', 'minor_AF', 'beta', 'se', 'pval', 'n_complete_samples']

    df2 = df2[df2['minor_AF'] >= 0.005]

    df2['pos'] = df2['pos'].apply(lambda x: int(x))
    df2['n_complete_samples'] = df2['n_complete_samples'].apply(lambda x: int(x))

    df2.to_csv(os.path.join(savedir, file), sep="\t", index=False)



def parallel_processing(sumstat_file_list, annotation_df, savedir, fnc, cores=1):
    pool = Pool(processes=cores)
    inputs = [(sumstat, annotation_df, savedir) for sumstat in sumstat_file_list]
    pool.starmap(fnc, inputs)
    pool.close()
    pool.join()


annotation_df = read_variant_annotation_file()
sumstat_file_list = []
for idx, row in df_metadata.iterrows():
    sumstat_path = os.path.join(ukb_gwas_dir, "{}.tsv".format(row['File']))
    sumstat_file_list.append(sumstat_path)

parallel_processing(sumstat_file_list, annotation_df, savedir, fnc=reformat_ukb_for_mr, cores=15)

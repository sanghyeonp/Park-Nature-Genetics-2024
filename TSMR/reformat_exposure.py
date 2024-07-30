import pandas as pd
import os
import code
# code.interact(local=dict(globals(), **locals()))

exposure_file = "MetS_noUKB_GWA.csv"

df = pd.read_csv(exposure_file, sep=",", index_col=False)

df = df[['SNP','CHR','BP','MAF','A1','A2','est','SE','Pval_Estimate']]
df['N'] = [1252787] * len(df)

df.to_csv("MetS_exposure.csv", sep=",", index=False)
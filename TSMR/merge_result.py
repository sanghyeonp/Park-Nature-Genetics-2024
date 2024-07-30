import pandas as pd
import os
import code
# code.interact(local=dict(globals(), **locals()))
import math
from math import isnan

result_dir = "./tsmr_result_individual"

results = [os.path.join(result_dir, f) for f in os.listdir(result_dir)]


df1 = None

for f in results:
    temp = pd.read_csv(f, sep=",", index_col=False)
    temp.replace("MR Egger", "MR-Egger", inplace=True)
    df1 = pd.concat([df1, temp])


df1.fillna("NA", inplace=True)

df1.to_csv("./MR_results_combined.tsv",
            sep="\t", index=False
            )
import os
import pandas as pd

files = [x for x in os.listdir(".") if x.startswith("query_chip_")]
dfc = pd.read_csv("query_chip_1.csv")
df1 = dfc["SNP"]
for file in files:
        dft = pd.read_csv(file)
        df1 = pd.merge(df1,dft, on='SNP', how='outer')
df1["Mean"] = df1.select_dtypes(include=['float64']).mean(axis=1).round(5)
df2 = df1.dropna()
df2 = df2.drop(df2[df2['Mean'] == 0].index)
df3 = df2[["SNP","Mean"]]
df2.to_csv("final_AC.csv",index=False)
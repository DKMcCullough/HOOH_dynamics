import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 


df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df = df_all
#df['log_sigma'] = df['log_sigma'].clip(lower=0) 
#df = df.rename(columns={"raw_abundance": "abundance"})

#inits = pd.read_csv("../data/inits/________.csv")

df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)] 
df_co = df.loc[df['assay'].str.contains('coculture', case=False)] 
df_mono = df.loc[~df['assay'].str.contains('coculture', case=False)] 
df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)] 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)] 


df = df_P

f,ax = plt.subplots(1,2,figsize=[8,5])

ax[0].plot(np.mean((df.rep1+df.rep3),axis=0),np.std((df.rep1+df.rep3),axis=0))
ax[1].plot(np.mean(np.((df.rep1+df.rep3),axis=0),np.std((df.rep1+df.rep3),axis=0))
plt.semilogy()
plt.show()



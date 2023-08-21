import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 


df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df = df_all


df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df.loc[df['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df.loc[~df['assay'].str.contains('coculture', case=False)].copy()  

df = df_mono 
df['techAmean'] = np.nanmean(np.r_[[df[i] for i in ['rep1','rep2']]],axis=0)
df['techBmean'] = np.nanmean(np.r_[[df[i] for i in ['rep3','rep4']]],axis=0)
df['techAstd'] = np.nanstd(np.r_[[df[i] for i in ['rep1','rep2']]],axis=0)
df['techBstd'] = np.nanstd(np.r_[[df[i] for i in ['rep3','rep4']]],axis=0)



#df['abundance_mean'] = np.nanmean(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0)


#df['log_abundance'] = np.nanmean(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))
#df['log_stdv'] = np.std(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))



df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 


df = df_P


f,ax = plt.subplots(1,2,figsize=[8,5])

ax[0].plot(np.mean((df.techAmean+df.techBmean),axis=0),np.std((df.techAmean+df.techBmean),axis=0))
ax[0].plot(np.mean((df.rep1+df.rep3)),axis=0),np.std(np.log((df.rep1+df.rep3),axis=0))
plt.semilogy()
plt.show()



#inits = pd.read_csv("../data/inits/________.csv")

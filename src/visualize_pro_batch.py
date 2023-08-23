import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 


df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df = df_all
'''

df['techAmean'] = np.nanmean(np.r_[[df[i] for i in ['rep1','rep2']]],axis=0)
df['techBmean'] = np.nanmean(np.r_[[df[i] for i in ['rep3','rep4']]],axis=0)
df['techAstd'] = np.nanstd(np.r_[[df[i] for i in ['rep1','rep2']]],axis=0)
df['techBstd'] = np.nanstd(np.r_[[df[i] for i in ['rep3','rep4']]],axis=0)
df['abundance_mean'] = np.nanmean(np.r_[[df[i] for i in ['techAmean','techBmean']]],axis=0)
df['abundance_std'] = np.nanstd(np.r_[[df[i] for i in ['techAmean','techBmean']]],axis=0)
#df_P['Pavg'] = df_P[['rep1', 'rep2', 'rep3', 'rep4']].mean(axis=1)
#df_P['Pstd'] = df_P[['rep1', 'rep2', 'rep3', 'rep4']].std(axis=1)

df['lr1'] = np.log(df['rep1'])
df['lr2'] = np.log(df['rep2'])
df['lr3'] = np.log(df['rep3'])
df['lr4'] = np.log(df['rep4'])


df['logtechA'] = np.nanmean(np.r_[[df[i] for i in ['lr1','lr2']]],axis=0)
df['logtechB'] = np.nanmean(np.r_[[df[i] for i in ['lr3','lr4']]],axis=0)
df['logtechAstd'] = np.nanstd(np.r_[[df[i] for i in ['lr1','lr2']]],axis=0)
df['logtechBstd'] = np.nanstd(np.r_[[df[i] for i in ['lr3','lr4']]],axis=0)
df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['logtechA','logtechB']]],axis=0)
df['log_abundance_std'] = np.nanstd(np.r_[[df[i] for i in ['logtechA','logtechB']]],axis=0)

'''
df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df.loc[df['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df.loc[~df['assay'].str.contains('coculture', case=False)].copy()  

#df = df_mono 

#df['log_abundance'] = np.nanmean(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))
#df['log_stdv'] = np.std(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))



df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 


df = df_P
#df =  df_S[df_S['strain'] == 'WH7803'].copy() 

'''

f1,ax = plt.subplots(1,2,figsize=[8,5],dpi = 300)
f1.suptitle('No logging  when Graphing')
ax[0].set_title('reps')
ax[1].set_title('total')
#ax[0].plot(np.log(df.techAmean),np.log(df.techAstd),label = 'a')
#ax[0].plot(np.log(df.techBmean),np.log(df.techBstd), label = 'b')
#ax[1].plot(np.log(df.abundance_mean),np.log(df.abundance_std), label = 'mean')
#ax[0].semilogy()
#ax[1].semilogy()
ax[0].plot((df.techAmean),(df.techAstd),label = 'a')
ax[0].plot((df.techBmean),(df.techBstd), label = 'b')
ax[1].plot((df.abundance_mean),(df.abundance_std), label = 'mean')
#ax[0].semilogy()
#ax[1].semilogy()

ax[0].legend()
ax[1].legend()

plt.show()

f1.savefig('../figures/f1',dpi=300)


f2,ax = plt.subplots(1,2,figsize=[8,5],dpi = 300)
f2.suptitle('Logged and then AVG')
ax[0].set_title('reps')
ax[1].set_title('total')
ax[0].plot((df.logtechA),(df.logtechAstd),label = 'a')
ax[0].plot((df.logtechB),(df.logtechBstd), label = 'b')
ax[1].plot((df.log_abundance),(df.log_abundance_std), label = 'mean')
#ax[0].semilogy()
#ax[1].semilogy()
ax[0].legend()
ax[1].legend()

plt.show()

f2.savefig('../figures/f2',dpi=300)
'''

strains = df_mono['strain'].unique()
nstrains = strains.shape[0]
colors = ['r','orange','g','b','k']
fig3, (ax)= plt.subplots(nstrains,2,figsize = (10,16))
fig3.suptitle('Unlogged Data')
for (S,si) in zip(strains,range(nstrains)):   
    df = df_mono[((df_mono['strain']==S))].copy()
    df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
    df['avg2'] = df[['rep2', 'rep4']].mean(axis=1)
    df['std1'] = df[['rep1', 'rep3']].std(axis=1)
    df['std2'] = df[['rep2', 'rep4']].std(axis=1)
    df0 = df[((df['assay']=='plus_0'))].copy()
    df400 = df[((df['assay']=='plus_400'))].copy()
    #print(df0.info(),df400.info())
    ax[si,0].plot(df0['time'],df0['avg1'], marker='o',label = S)
    ax[si,1].plot(df400['time'],df400['avg1'], yerr='std1',marker='o',label = S)
    df0.plot(kind='scatter', x='time', y ='avg1', yerr='std1',style="-", label = S, title = '0 HOOH assay', ylabel = 'cells per mL',logy = True)
    df400.plot(kind='scatter', x='time', y ='avg1', yerr='std1',style="-", label = S, title = '400 HOOH assay', ylabel = 'cells per mL',logy = True)



#inits = pd.read_csv("../data/inits/________.csv")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 


df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df = df_all

df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

#df = df_mono 

#df['log_abundance'] = np.nanmean(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))
#df['log_stdv'] = np.std(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))



df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 


df = df_mono

plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'

strains = df_mono['strain'].unique()
nstrains = strains.shape[0]
colors = ('orange', 'r', 'green', 'c', 'purple', 'k')
#need to slice by Vol number !!! (2 cat + Syns)
fig1, (ax0,ax1)= plt.subplots(1,2,figsize = (16,9))
fig1.suptitle('Pro and Syn Monocultures')

for (S,si) in zip(strains,range(nstrains)):   
    count = si
    df = df_mono[((df_mono['strain']==S))].copy()
    df['log1'] = np.log(df['rep1'])
    df['log2'] = np.log(df['rep2'])
    df['log3'] = np.log(df['rep3'])
    df['log4'] = np.log(df['rep4'])
    df['avg1'] = df[['log1', 'log3']].mean(axis=1)
    df['avg2'] = df[['log2', 'log4']].mean(axis=1)
    df['std1'] = df[['log1', 'log3']].std(axis=1)
    df['std2'] = df[['log2', 'log4']].std(axis=1)
    df0 = df[((df['assay']=='plus_0'))].copy()
    df400 = df[((df['assay']=='plus_400'))].copy()
    ax0.errorbar(df0['time'],df0['avg1'],yerr=df0['std1'], marker='o',color = colors[count], label = str(S)+' avg 1')
    ax0.errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='v',color = colors[count], label = str(S)+' avg 2')
    ax1.errorbar(df400['time'],df400['avg1'], yerr=df400['std1'],marker='o', color = colors[count], label = str(S)+' avg 1')
    ax1.errorbar(df400['time'],df400['avg2'], yerr=df400['std2'],marker='v', color = colors[count], label = str(S)+' avg 2')
    ax1.set_ylim(6.5,15)
    ax0.set_ylim(6.5,15)
    l1  = ax0.legend(loc = 'lower left')
    l1.draw_frame(False)
    
# make space on the right for annotation (e.g. ROS=0, etc.)
fig1.subplots_adjust(right=0.85, wspace = 0.25, hspace = 0.25)

# titles and labels 
ax0.set_title('Monocultures in 0 HOOH')
ax1.set_title('Monocultures in 400 HOOH')

ax0.set_xlabel('Time (days)')
ax1.set_xlabel('Time (days)')
ax0.set_ylabel('Cells (ml$^{-1}$)')
ax1.set_ylabel('Cells (ml$^{-1}$)')

fig1.savefig('../figures/mono_all_graphed')
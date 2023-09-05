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


df = df_P
#df =  df_S[df_S['strain'] == 'WH7803'].copy() 


#####################################################
# graphing monoculture data by strain in unlogged df
#####################################################

strains = df_mono['Vol_number'].unique()
nstrains = strains.shape[0]
#need to slice by Vol number !!! (2 cat + Syns)
fig3, (ax)= plt.subplots(nstrains,2,figsize = (12,14))
fig3.suptitle('Unlogged Data', size =25 )
for (S,si) in zip(strains,range(nstrains)):   
    df = df_mono[((df_mono['Vol_number']==S))].copy() #setting df to only monocultre data that share the vol number.
    df['avg1'] = df[['rep1', 'rep3']].mean(axis=1) #avging over odd or even reps for bio reps of assay s (1/2 and 3/4 are technical reps) 
    df['avg2'] = df[['rep2', 'rep4']].mean(axis=1)
    df['std1'] = df[['rep1', 'rep3']].std(axis=1)
    df['std2'] = df[['rep2', 'rep4']].std(axis=1)
    df0 = df[((df['assay']=='plus_0'))].copy() #split df inot 0 and 400 assays 
    df400 = df[((df['assay']=='plus_400'))].copy()
    #print(df0.info(),df400.info())
    ax[si,0].errorbar(df0['time'],df0['avg1'],yerr=df0['std1'], marker='o',label = 'vol num ' + str(S)+'  avg 1')
    ax[si,0].errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='v',label = 'vol num ' + str(S)+'  avg 2')
    ax[si,1].errorbar(df400['time'],df400['avg1'], yerr=df400['std1'],marker='o',label ='vol num ' + str(S)+'  avg 1')
    ax[si,1].errorbar(df400['time'],df400['avg2'],yerr=df400['std2'], marker='v',label = 'vol num ' + str(S)+'  avg 2')
    #ax[si,0].set_ybound(1000,2500000)
    #ax[si,1].set_ybound(1000,1500000)
    l3  = ax[si,1].legend(loc = 'upper center')
    l3.draw_frame(False)
    #df0.plot(kind='scatter', x='time', y ='avg1', yerr='std1',style="-", label = S, title = '0 HOOH assay', ylabel = 'cells per mL',logy = True)
    #df400.plot(kind='scatter', x='time', y ='avg1', yerr='std1',style="-", label = S, title = '400 HOOH assay', ylabel = 'cells per mL',logy = True)


# make space on the right for annotation (e.g. ROS=0, etc.)
fig3.subplots_adjust(right=0.85, wspace = 0.25, hspace = 0.25)

# titles
ax[0,0].set_title('Monocultures in 0 HOOH')
ax[0,1].set_title('Monocultures in 400 HOOH')

# xlabels
for a in ax[-1,:]:
    a.set_xlabel('Time (days)')

# ylabels
for a in ax[:,0]:
    a.set_ylabel('Cells (ml$^{-1}$)')

fig3.savefig('../figures/monoculture_graphs')
#inits = pd.read_csv("../data/inits/________.csv")


#####################################################
# graphing logged data 
#####################################################



fig4, (ax)= plt.subplots(nstrains,2,figsize = (12,14))
fig4.suptitle('Logged Data',size =25)
for (S,si) in zip(strains,range(nstrains)):    #zipping the strain list with number of strains to have loop for whole df
    df = df_mono[((df_mono['Vol_number']==S))].copy() #setting df to that which matches vol number in Strain list at said itteration 
    df['log1'] = np.log(df['rep1']) #logging and stdv for data and error evalution 
    df['log2'] = np.log(df['rep2']) #logging reps
    df['log3'] = np.log(df['rep3'])
    df['log4'] = np.log(df['rep4'])
    df['avg1'] = df[['log1', 'log3']].mean(axis=1) #avg of logged data
    df['avg2'] = df[['log2', 'log4']].mean(axis=1) 
    df['std1'] = df[['log1', 'log3']].std(axis=1)
    df['std2'] = df[['log2', 'log4']].std(axis=1)
    df0 = df[((df['assay']=='plus_0'))].copy() #seleection only assay 0 from working df
    df400 = df[((df['assay']=='plus_400'))].copy() #selecting 400 H assay from loop's working df
    ax[si,0].errorbar(df0['time'],df0['avg1'],yerr=df0['std1'], marker='o',label = 'vol num ' + str(S)+'  avg 1')
    ax[si,0].errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='v',label ='vol num ' + str(S)+'  avg 2')
    ax[si,1].errorbar(df400['time'],df400['avg1'], yerr=df400['std1'],marker='o',label = 'vol num ' + str(S)+' avg 1')
    ax[si,1].errorbar(df400['time'],df400['avg2'],yerr=df400['std2'], marker='v',label = 'vol num ' + str(S)+' avg 2')
    ax[si,0].semilogy()
    ax[si,1].semilogy()
    l4  = ax[si,1].legend(loc = 'upper center')
    l4.draw_frame(False)

#Config fig
# make space on the right for annotation (e.g. ROS=0, etc.)
fig4.subplots_adjust(right=0.85, wspace = 0.25, hspace = 0.25)

# titles
ax[0,0].set_title('Monocultures in 0 HOOH')
ax[0,1].set_title('Monocultures in 400 HOOH')

# xlabels
for a in ax[-1,:]:
    a.set_xlabel('Time (days)')

# ylabels
for a in ax[:,0]:
    a.set_ylabel('Cells (ml$^{-1}$)')

plt.show() 

#fig4.savefig('../figures/monoculture_logged_graphs')
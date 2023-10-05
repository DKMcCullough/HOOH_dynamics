'''

name:   model_abiotic_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model of Monoculture BCC assays to graph 0 H phyotplankton biomass 

working on: - getting this in model to play so we can model all at once 

'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 


#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'

######################################################
#reading in data and configureing 
#####################################################

df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df = df_all

#make unlogged avgs
df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
df['std1'] = df[['rep1', 'rep3']].std(axis=1)

#make logs of data
df['log1'] = np.log(df['rep1'])
df['log3'] = np.log(df['rep3'])

#avg and std of logged data
df['lavg1'] = df[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
df['stdlog1'] = df[['log1', 'log3']].std(axis=1) #taking stdv of logged reps

#splicing abiotic and mono or coculture data
df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df.loc[df['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df.loc[~df['assay'].str.contains('coculture', case=False)].copy()  
df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 

#setting working df
df = df_abiotic
#picking 0 or 400 assay dfs
df0 = df.loc[df['assay'].str.contains('_0', case=False)].copy()
df400 = df.loc[df['assay'].str.contains('_0', case=False)].copy()

#####################################################
# graphing abiotic data by assay in unlogged df
#####################################################

fig1, (ax1)= plt.subplots(2,1,sharex=True,figsize = (10,8))
fig1.suptitle('Unlogged Abiotic Data', size =25 )
fig1.supylabel('HOOH concentration (nM)')
fig1.supxlabel('Time (Days)')
ax1[0].set_title('Monocultures in 0 HOOH')
ax1[1].set_title('Monocultures in 400 HOOH')

ax1[0].semilogy()
ax1[1].semilogy()

# make space on the right for annotation (e.g. ROS=0, etc.)
fig1.subplots_adjust(right=0.90, left=0.20,wspace = 0.25, hspace = 0.30)


#graphing data
ax1[0].errorbar(df0['time'],df0['avg1'],yerr=df0['std1'], marker='o',label = 'avg unlogged')
ax1[0].plot(df0['time'],df0['rep1'], label = 'rep1')
ax1[0].plot(df0['time'],df0['rep3'], label = 'rep3')
ax1[1].errorbar(df400['time'],df400['avg1'], yerr=df400['std1'],marker='o',label ='avg unlogged')
ax1[1].plot(df400['time'],df400['rep1'], label = ' rep1')
ax1[1].plot(df400['time'],df400['rep3'], label = 'rep3')
#config legend 
l1  = ax1[1].legend(loc = 'upper left')
l1.draw_frame(False)

plt.show()
#save fig 
fig1.savefig('../figures/abiotic_unlogged')



#####################################################
# graphing abiotic data by assay in logged df
#####################################################

fig2, (ax2)= plt.subplots(2,1,sharex=True,figsize = (10,8))
fig2.suptitle('Logged Abiotic Data', size =25 )
fig2.supylabel('HOOH concentration (nM)')
fig2.supxlabel('Time (Days)')
ax1[0].set_title('Abitoic 0 HOOH Spike Data')
ax2[1].set_title('Abitoic 400 HOOH Spike Data')

ax2[0].semilogy()
ax2[1].semilogy()

# make space on the right for annotation (e.g. ROS=0, etc.)
fig2.subplots_adjust(right=0.90, left=0.20, wspace = 0.25, hspace = 0.35)


#graphing data
ax2[0].errorbar(df0['time'],df0['lavg1'],yerr=df0['stdlog1'], marker='o',label = 'avg loggged')
ax2[0].plot(df0['time'],df0['log1'], label = 'log1')
ax2[0].plot(df0['time'],df0['log3'], label = 'log3')
ax2[1].errorbar(df400['time'],df400['lavg1'], yerr=df400['stdlog1'],marker='o',label ='avg logged')
ax2[1].plot(df400['time'],df400['log1'], label = 'log1')
ax2[1].plot(df400['time'],df400['log3'], label = 'log3')
#config legend 
l2  = ax2[1].legend(loc = 'upper left')
l2.draw_frame(False)

plt.show()

#save fig 

fig2.savefig('../figures/abiotic_logged')

###############################
#graph Stats of 0 assay  
######################
fig3,ax3 = plt.subplots(1,2,figsize=[10,8],dpi = 300)
fig3.suptitle('(0 H spike Assay )')
#config 
ax3[0].set_xlabel('unlogged mean')
ax3[0].set_ylabel('unllogged stdv')
ax3[1].set_xlabel('logged mean')
ax3[1].set_ylabel('standard deviation')
ax3[0].semilogy()
ax3[1].semilogy()

fig3.subplots_adjust(right=0.90, left=0.20, wspace = 0.25, hspace = 0.35)

#plot data

ax3[0].plot(df0.avg1,df0.std1,label = 'unlogged stats')
ax3[1].plot(df0.lavg1,df0.stdlog1,label = 'Logged stats')

l3  = ax3[1].legend(loc = 'upper left')
l3.draw_frame(False)
l0  = ax3[0].legend(loc = 'upper left')
l0.draw_frame(False)

plt.show()

fig3.savefig('../figures/0_abiotic_stats',dpi=300)

#Stats of 400 assay  

fig4,ax4 = plt.subplots(1,2,figsize=[10,8],dpi = 300)
fig4.suptitle('400 H spike Assay ')

ax4[0].set_xlabel('unlogged mean')
ax4[0].set_ylabel('unlogged stdv')
ax4[1].set_xlabel('mean')
ax4[1].set_ylabel('standard deviation')
ax4[0].semilogy()
ax4[1].semilogy()


ax4[0].plot(df400.avg1,df400.std1,label = 'unlogged stats')
ax4[1].plot(df400.lavg1,df400.stdlog1,label = 'Logged stats')

l4  = ax4[1].legend(loc = 'upper left')
l4.draw_frame(False)
l5  = ax4[0].legend(loc = 'upper left')
l5.draw_frame(False)

plt.show()

fig4.savefig('../figures/400_abiotic_stats',dpi=300)






fig4.savefig('../figures/400_abiotic_stats',dpi=300)



















print("I'm done with Abiotic Spike Assays bro! ")
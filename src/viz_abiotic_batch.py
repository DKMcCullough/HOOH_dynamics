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

df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['log4'] = np.log(df['rep4'])

#####

#bio rep avgs put together 
df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
df['avg2'] = df[['rep2', 'rep4']].mean(axis=1)
df['std1'] = df[['rep1', 'rep3']].std(axis=1)
df['std2'] = df[['rep2', 'rep4']].std(axis=1)
#bio rep logged avgs and stdvs 
df['lavg1'] = df[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
df['lavg2'] = df[['log2', 'log4']].mean(axis=1)
df['stdlog1'] = df[['log1', 'log3']].std(axis=1) #taking stdv of logged reps
df['stdlog2'] = df[['log2', 'log4']].std(axis=1)


#total avgs and stdvs
df['abundance'] =  np.nanmean(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0)
df['sigma'] = np.std(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0)
df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['log1','log2','log3','log4']]],axis=0)
df['log_sigma'] = np.std(np.r_[[df[i] for i in ['log1','log2','log3','log4']]],axis=0)


#splicing abiotic and mono or coculture data
df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df.loc[df['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df.loc[~df['assay'].str.contains('coculture', case=False)].copy()  
df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[~df_mono['organism'].str.contains('P', case=False)].copy() 

#setting working df
df = df_abiotic
#picking 0 or 400 assay dfs
df0 = df.loc[df['assay'].str.contains('_0', case=False)].copy()
df400 = df.loc[df['assay'].str.contains('_0', case=False)].copy()

treats = df['assay'].unique()
ntreats = treats.shape[0]


    #setting up fig 1 for dyanmics 


#####################################################
# Set up large loop  (treatments) 
#####################################################

for (t,nt) in zip(treats,range(ntreats)):

    df = df_abiotic[((df_abiotic['assay']==t))].copy()
    
    fig1,ax1 = plt.subplots(1,2,figsize=[11,8]) #plot creation and config 
    fig1.suptitle('Raw dynamics of'+ str(t)) #full title config
    ax1[0].set_title('HOOH Concentration')
    ax1[1].set_title('avg vs stdv')
    fig1.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.35, hspace=0.30) #shift white space for better fig view
    ax1[0].set_xlabel('Time (days)')
    ax1[0].set_ylabel('HOOH concentration (\u03BCM)')
    ax1[1].set_xlabel('Raw Mean')
    ax1[1].set_ylabel('Raw standard deviation ')

    fig2,ax2 = plt.subplots(1,2,figsize=[11,8]) #plot creation and config 
     #full title config
    fig2.suptitle('Logged dynamics of '+ str(t))
    ax2[0].set_title('HOOH concentration ')
    ax2[1].set_title('Log avg vs stdv')
    fig2.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.30, hspace=0.30) #shift white space for better fig view
    ax2[0].set_xlabel('Time (days)')
    ax2[0].set_ylabel('HOOH concentration (\u03BCM)')
    ax2[1].set_xlabel('Log Mean')
    ax2[1].set_ylabel('Log Standard deviation')


    
#####################################################
# Set up loop of vol numbers inside treatment loop  
#####################################################


#setting working df as a single Experiment in df_all
    ax1[0].plot(df.time,df.rep1, label = 'rep1')
    ax1[0].plot(df.time,df.rep2, label = 'rep2')
    ax1[0].plot(df.time,df.rep3, label = 'rep3')
    ax1[0].plot(df.time,df.rep4, label = 'rep4')
    ax1[0].errorbar(df.time,df.abundance, yerr=df.sigma, marker = '*', c='g',label =  'Mean')
    ax1[0].errorbar(df.time,df.avg1, yerr=df.std1, marker = 'o', c='b',label =  'Avg 1')
    ax1[0].errorbar(df.time,df.avg2, yerr=df.std2, marker = 'd', c='r',label =  'Avg 2')
    ax1[1].scatter(df.abundance,df.sigma, c='g')
    ax1[1].scatter(df.avg1,df.std1, c='b')
    ax1[1].scatter(df.avg2,df.std2, c='r')
    ax1[0].semilogy()
    l1 = ax1[0].legend(loc = 'upper left')
    l1.draw_frame(False)
    ax2[0].plot(df.time,df.log1, label = 'l1')
    ax2[0].plot(df.time,df.log2, label = 'l2')
    ax2[0].plot(df.time,df.log3, label = 'l3')
    ax2[0].plot(df.time,df.log4, label = 'l4')
    ax2[0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma, marker = '*', c='g',label =  'Log Mean')
    ax2[0].errorbar(df.time,df.lavg1, yerr=df.stdlog1, marker = 'o', c='b',label =  'Log Avg 1')
    ax2[0].errorbar(df.time,df.lavg2, yerr=df.stdlog2, marker = 'd', c='r',label =  'Log Avg 2')
    ax2[1].scatter(df.log_abundance,df.log_sigma, c='g')
    ax2[1].scatter(df.lavg1,df.stdlog1, c='b')
    ax2[1].scatter(df.lavg2,df.stdlog2, c='r')
    ax2[0].semilogy()
    ax2[1].semilogy()
    ax2[1].semilogx()
    l2 = ax2[0].legend(loc = 'upper left')
    l2.draw_frame(False)








plt.show











print("I'm done with Abiotic Spike Assays bro! ")
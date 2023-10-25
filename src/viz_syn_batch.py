'''

name:   vizualize_syn_batch.py 

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
df_all.fillna(0)


df = df_all


#logging data for latter graphing 
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

#slicing into different treatments and organisms 
df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 

#####################################################
# graphing monoculture data by strain in unlogged df
#####################################################

df = df_S
strains = df['Vol_number'].unique()
strains = strains[~pd.isna(strains)] #getting rid og annoying nan in list 

nstrains = strains.shape[0] #making numnber list for exp

treats = df['assay'].unique()
ntreats = treats.shape[0]




#####################################################
# Set up large loop  (treatments) 
#####################################################

for (t,nt) in zip(treats,range(ntreats)):

    df_w = df_S[((df_S['assay']==t))].copy()
    #setting up fig 1 for dyanmics 
    fig1,ax1 = plt.subplots(nstrains,2,figsize=[12,12]) #plot creation and config 
    fig1.suptitle('Syn Data in '+ str(t)) #full title config
    ax1[0,0].set_title('Syn Cell Concentration')
    fig1.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.30, hspace=0.30) #shift white space for better fig view


    fig2,ax2 = plt.subplots(nstrains,2,figsize=[12,12]) #plot creation and config 
    fig2.suptitle('Syn Logged Data in '+ str(t)) #full title config
    ax2[0,0].set_title('Logged Syn Cell Concentration')

    fig2.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.30, hspace=0.30) #shift white space for better fig view

#####################################################
# Set up loop of vol numbers inside treatment loop  
#####################################################

    for (S,si) in zip(strains,range(nstrains)):   
        df = df_w[((df_w['Vol_number']==S))].copy()
        #setting working df as a single Experiment in df_all
        ax1[si,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker = '*', c='orange',label =  'Avg ' +str(S) ) 
        ax1[si,0].errorbar(df.time,df.avg1, yerr=df.std1, marker = 'o', c='b',label =  'Avg 1 ' ) 
        ax1[si,0].errorbar(df.time,df.avg2, yerr=df.std2, marker = 'o', c='r',label =  'Avg 2 ') 
        ax1[si,1].scatter(df.abundance,df.sigma, c='k')
        ax1[si,1].scatter(df.avg1,df.std1, c='b')
        ax1[si,1].scatter(df.avg2,df.std2, c='r')
        ax1[si,0].semilogy()
        l1 = ax1[si,0].legend(loc = 'upper left')
        l1.draw_frame(False)
        ax2[si,0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma,  marker = '*', c='darkred',label =  'Log Avg ' +str(S) ) 
        ax2[si,0].errorbar(df.time,df.lavg1, yerr=df.stdlog1,  marker = 'o', c='k',label =  'Log Avg 1 ' ) 
        ax2[si,0].errorbar(df.time,df.lavg2, yerr=df.stdlog2, marker = 'o', c='r',label =  'Log Avg 2 ') 
        ax2[si,1].scatter(df.log_abundance,df.log_sigma, c='k')
        ax2[si,1].scatter(df.lavg1,df.stdlog1, c='b')
        ax2[si,1].scatter(df.lavg2,df.stdlog2, c='r')
        ax2[si,0].semilogy()
        ax2[si,1].semilogy()
        ax2[si,1].semilogx()
        l2 = ax2[si,0].legend(loc = 'upper left')
        l2.draw_frame(False)
        
        '''
        #ax1 labels
    # xlabels
        for a in ax1[-1,1]:
            a.set_xlabel('Time (days)')

    # ylabels
        for a in ax1[:,0]:
            a.set_ylabel('Cells (ml$^{-1}$)')
        
        for a in ax1[-1,1]:
            a.set_ylabel('STDV')
            a.set_xlabel('Abundance')
        
        
        
        #ax2 labels 
        for a in ax2[-1,1]:
            a.set_xlabel('Time (days)')

    # ylabels
        for a in ax2[:,0]:
            a.set_ylabel('Cells (ml$^{-1}$)')
            
        for a in ax2[-1,1]:
            a.set_ylabel('STDV')
            a.set_xlabel('Abundance')
            '''

plt.show() 
#fig1.savefig('../figures/Syn_raw_graphs')
#fig2.savefig('../figures/Syn_logged_graphs')

# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print(' Done my guy ')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



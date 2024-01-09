'''

name:   viz_pro_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model of Monoculture BCC assays to graph Pro  phyotplankton biomass 



'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 



######################################################
#reading in data and configureing 
#####################################################

#df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all = pd.read_csv("../data/BCC_2-5-dataset.csv",header=1)
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
df_S = df_mono.loc[~df_mono['organism'].str.contains('P', case=False)].copy() 


df = df_P


treats = df['assay'].unique()
ntreats = treats.shape[0]


    #setting up fig 1 for dyanmics 


#####################################################
# Set up large loop  (treatments) 
#####################################################

for (t,nt) in zip(treats,range(ntreats)):

    df = df_P[((df_P['assay']==t))].copy()
    
    fig1,ax1 = plt.subplots(1,2,figsize=[11,8]) #plot creation and config 
    fig1.suptitle('Pro Data in '+ str(t)) #full title config
    ax1[0].set_title('Pro Cell Concentration')
    ax1[1].set_title('avg vs stdv')
    fig1.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.35, hspace=0.30) #shift white space for better fig view


    fig2,ax2 = plt.subplots(1,2,figsize=[11,8]) #plot creation and config 
     #full title config
    fig2.suptitle('Pro Logged Data in '+ str(t))
    ax2[0].set_title('Pro Cell Concentration')
    ax2[1].set_title('avg vs stdv')
    fig2.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.30, hspace=0.30) #shift white space for better fig view

    ax1[0].set_xlabel('Time (days)')
    ax1[0].set_ylabel('Cell Concentration mL-1')
    ax1[1].set_xlabel('standard deviation')
    ax1[1].set_ylabel('mean value')
    ax2[0].set_xlabel('Time (days)')
    ax2[0].set_ylabel('Cell Concentration mL-1')
    ax2[1].set_xlabel('Average')
    ax2[1].set_ylabel('Standard deviation')
    
    fig3,ax3 = plt.subplots(1,2,figsize = [8,6])
    fig3.suptitle("Pearsons's R coorelations "+ str(t))
    fig3.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.85, wspace=0.25, hspace=0.4)
    ax3[0].set_title('Raw ')
    ax3[1].set_title('Logged ')
#####################################################
# Set up loop of vol numbers inside treatment loop  
#####################################################


#setting working df as a single Experiment in df_all
    ax1[0].errorbar(df.time,df.abundance, yerr=df.sigma, marker = '*', c='g',label =  'Mean')
    ax1[0].errorbar(df.time,df.avg1, yerr=df.std1, marker = 'o', c='b',label =  'Avg 1')
    ax1[0].errorbar(df.time,df.avg2, yerr=df.std2, marker = 'd', c='r',label =  'Avg 2')
    ax1[1].scatter(df.abundance,df.sigma, c='g')
    ax1[1].scatter(df.avg1,df.std1, c='b')
    ax1[1].scatter(df.avg2,df.std2, c='r')
    ax1[0].semilogy()
    l1 = ax1[0].legend(loc = 'center')
    l1.draw_frame(False)
    ax2[0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma, marker = '*', c='g',label =  'Log Mean')
    ax2[0].errorbar(df.time,df.lavg1, yerr=df.stdlog1, marker = 'o', c='b',label =  'Log Avg 1')
    ax2[0].errorbar(df.time,df.lavg2, yerr=df.stdlog2, marker = 'd', c='r',label =  'Log Avg 2')
    ax2[1].scatter(df.log_abundance,df.log_sigma, c='g')
    ax2[1].scatter(df.lavg1,df.stdlog1, c='b')
    ax2[1].scatter(df.lavg2,df.stdlog2, c='r')
    ax2[0].semilogy()
    ax2[1].semilogy()
    ax2[1].semilogx()
    l2 = ax2[0].legend(loc = 'center')
    l2.draw_frame(False)

    raw_r,raw_p  = scipy.stats.pearsonr(df['abundance'],df['sigma'])
    log_r,log_p = scipy.stats.pearsonr(df['log_abundance'],df['log_sigma'])
    #print(raw_r,log_r)
    ax3[0].hist(raw_r,color = 'red')
    ax3[1].hist(log_r,color = 'b')


fig1.savefig('../figures/Pro_'+(str(t))+'raw_graphs')
fig2.savefig('../figures/Pro_'+(str(t))+'logged_graphs')
fig3.savefig('../figures/Pro_'+(str(t))+'corr_graphs')


plt.show() 


# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print(' Done my guy ')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')




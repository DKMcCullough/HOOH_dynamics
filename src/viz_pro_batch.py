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














'''
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
    ax[si,0].semilogy()
    ax[si,1].semilogy()
    l3  = ax[si,1].legend(loc = 'upper center')
    l3.draw_frame(False)


# make space on the right for annotation (e.g. ROS=0, etc.)
fig3.subplots_adjust(right=0.90,wspace = 0.15, hspace = 0.25)

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
fig4.subplots_adjust(right=0.90,wspace = 0.45, hspace = 0.25)

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

fig4.savefig('../figures/monoculture_logged_graphs')
'''



plt.show() 


# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print(' Done my guy ')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')




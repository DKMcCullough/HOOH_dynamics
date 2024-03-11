'''

name:   viz_buffer.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser/Projects/ROS.../HOOH.....src'

author: DKM



'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import seaborn as sns

#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'

######################################################
#reading in data and configureing 
#####################################################
df_spike = pd.read_excel("../data/ROS_uncertainty.xlsx",sheet_name = 'Spike_assays', header = 0)
df_media = pd.read_excel("../data/ROS_uncertainty.xlsx",sheet_name = 'Media_assays', header = 0)
df_growth = pd.read_excel("../data/ROS_uncertainty.xlsx",sheet_name = 'Growth_assays', header = 0)
df_buffer = pd.read_excel("../data/ROS_uncertainty.xlsx",sheet_name = 'Buffer_assays', header = 0)

df_all =   df_buffer
df = df_all

df = df.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df.fillna(0)
df[['A_tech_rep1','A_tech_rep2','B_tech_rep1','B_tech_rep2','C_tech_rep1','C_tech_rep2']].astype('float64')


#logging data for latter graphing 
df['logA1'] = np.log(df['A_tech_rep1'])
df['logA2'] = np.log(df['A_tech_rep2'])
df['logB1'] = np.log(df['B_tech_rep1'])
df['logB2'] = np.log(df['B_tech_rep2'])
df['logC1'] = np.log(df['C_tech_rep1'])
df['logC2'] = np.log(df['C_tech_rep2'])


#####

#bio rep avgs put together 
df['avgA'] = df[['A_tech_rep1', 'A_tech_rep2']].mean(axis=1)
df['avgB'] = df[['B_tech_rep1', 'B_tech_rep2']].mean(axis=1)
df['avgC'] = df[['C_tech_rep1', 'C_tech_rep2']].mean(axis=1)

df['stdA'] = df[['A_tech_rep1', 'A_tech_rep2']].std(axis=1)
df['stdB'] = df[['B_tech_rep1', 'B_tech_rep2']].std(axis=1)
df['stdC'] = df[['C_tech_rep1', 'C_tech_rep2']].std(axis=1)


#bio rep logged avgs and stdvs 
df['lavgA'] = df[['logA1', 'logA2']].mean(axis=1)
df['lavgB'] = df[['logB1', 'logB2']].mean(axis=1)
df['lavgC'] = df[['logC1', 'logC2']].mean(axis=1)

df['stdlogA'] = df[['logA1', 'logA2']].std(axis=1) #taking stdv of logged reps
df['stdlogB'] = df[['logB1', 'logB2']].std(axis=1)
df['stdlogC'] = df[['logC1', 'logC2']].std(axis=1)


#total avgs and stdvs
df['abundance'] =  np.nanmean(np.r_[[df[i] for i in ['A_tech_rep1', 'A_tech_rep2','B_tech_rep1', 'B_tech_rep2','C_tech_rep1', 'C_tech_rep2']]],axis=0)
df['sigma'] = np.nanstd(np.r_[[df[i] for i in ['A_tech_rep1', 'A_tech_rep2','B_tech_rep1', 'B_tech_rep2','C_tech_rep1', 'C_tech_rep2']]],axis=0)
df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['logA1','logA2','logB1','logB2','logC1','logC2']]],axis=0)
df['log_sigma'] = np.nanstd(np.r_[[df[i] for i in ['logA1','logA2','logB1','logB2','logC1','logC2']]],axis=0)


#graph all at once 

#sns.relplot(data=df, x="abundance", y="sigma", palette = 'cool',size ='organism', hue="Buffer",style= 'Light', markers =True,  kind="scatter").set(title='Raw Data' )
sns.relplot(data=df, x="abundance", y="sigma", palette = 'cool',size ='Buffer_Concentration (microM) ', hue="Buffer",style= 'Light', markers =True,  kind="scatter").set(title='Raw Data' )

#sns.relplot(data=df, x="log_abundance", y="log_sigma", palette = 'cool',size ='organism', hue="Buffer", style= 'Light', markers =True,  kind="scatter").set(title='Log Data' )
sns.relplot(data=df, x="log_abundance", y="log_sigma", palette = 'cool',size ='Buffer_Concentration (microM) ', hue="Buffer",style= 'Light', markers =True,  kind="scatter").set(title='Log Data' )



#make arrays for going through loop via assay 

orgs = df['organism'].unique()
norgs = orgs.shape[0]

dfs = df['Buffer'].unique()
ndfs = dfs.shape[0]


#####################################################
# Set up large loop  for graphing 
#####################################################

for d,nd in zip(orgs,range(norgs)):
    dfw = df[(df['organism'] == d) ]



    a = sns.relplot(data=dfw, x="abundance", y="sigma", palette = 'cool',size ='Buffer_Concentration (microM) ', hue='Buffer', style = 'Light', markers =True,  kind="scatter").set(title='Raw '+str(d)+' Data' )
    
    a.ax.set_xlabel("Raw Mean HOOH",fontsize=13)
    a.ax.set_ylabel("Raw Sigma",fontsize=13)
    a.ax.tick_params(labelsize=12)
    
    b = sns.relplot(data=dfw, x="log_abundance", y="log_sigma", palette = 'cool',size ='Buffer_Concentration (microM) ', hue='Buffer', style = 'Light', markers =True,  kind="scatter").set(title='Log '+str(d)+' Data')
    b.ax.set_xlabel("Log Mean HOOH",fontsize=13)
    b.ax.set_ylabel("Log Sigma",fontsize=13)
    b.ax.tick_params(labelsize=12)



    plt.show() 
    
    a.savefig('../figures/raw_'+str(d)+'data')
    b.savefig('../figures/Log_'+str(d)+'data')


# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print(' Done my guy ')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')




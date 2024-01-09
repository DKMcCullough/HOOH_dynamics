
'''

name:   model_pro_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/HOOH_dynamics/src'
    
author: DKM

goal: Loop model and graph 0 H Pro assay and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop for 0 and 400 using different init files? (need H connected to 0 H first) 

'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ODElib
import random as rd
import sys


######################################################
#reading in data and configureing 
#####################################################
df_1 = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_1-31-dataset', header = 1)
df_2 = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_2-5-dataset', header = 1)

df_all = df_1


#df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#slicing data into abiotic, biotic, and Pro only dataframes

df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 

#setting working df as pro only 
df = df_P

#####################################################
#config data in df from raw for odelib usefulness
#####################################################

#making avg columns of technical reps (std here only for graphing, not logged here)
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['log4'] = np.log(df['rep4'])
df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
df['avg2'] = df[['rep2', 'rep4']].mean(axis=1)
df['std1'] = df[['rep1', 'rep3']].std(axis=1)
df['std2'] = df[['rep2', 'rep4']].std(axis=1)

df['lavg1'] = df[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
df['lavg2'] = df[['log2', 'log4']].mean(axis=1)
df['stdlog1'] = df[['log1', 'log3']].std(axis=1) #taking stdv of logged reps
df['stdlog2'] = df[['log2', 'log4']].std(axis=1)

#setting working df for model as far as abundance and log abundance values 
df.rename(columns = {'avg1':'abundance'}, inplace = True) #reaneme main df column to be fit by odelib 
df.rename(columns = {'lavg1':'log_abundance'}, inplace = True) #reaneme log of main df column to be fit by odelib 

#splitting df of Pro into 0 and 400 H assays 
df0 = df.loc[~ df['assay'].str.contains('4', case=False)]  #assay 0 H 
df4 = df.loc[df['assay'].str.contains('4', case=False)]  #assay 400 H (actually around 360 nM in data)

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig2.suptitle('Pro  Monocultures')

#format fig  
ax0.set_title('Pro in 0 HOOH ') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('Pro in 400 HOOH ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)')#settign x axis label for graph 2 

#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df0['time'],df0['abundance'],yerr=df0['std1'], marker='o', label = 'avg1')
ax0.errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='v', label = 'avg2')
# graph 400 H assay even and odd avgs
ax1.errorbar(df4['time'],df4['abundance'],yerr=df4['std1'], marker='o', label = 'avg1')
ax1.errorbar(df4['time'],df4['avg2'],yerr=df4['std2'], marker='v', label = 'avg2')

plt.legend()

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits0 = pd.read_csv("../data/inits/pro9215_inits0.csv")

#setting how many MCMC chains you will run 
nits = 10000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['P','N'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.000002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
#setting state variiable  prior guess
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+5})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+6})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
P0_mean = 50000
N0_mean = 900000

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_0H,
                          parameter_names=['k1','k2','P0','N0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          P0 = P0_prior.copy(),
                          N0  = N0_prior.copy(),
                          t_steps=1000,
                          P = P0_mean,
                          N = N0_mean,
                            )
    return M

def mono_0H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2 = params[0], params[1]
    P,N = y[0],y[1]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P     
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    return [dPdt,dNdt]

#df0.loc[:,'log_abundance'] = np.log(10**df0.log_abundance)

# get_models
a0 = get_model(df0) 

#broken here!!!!!!!!!!
# do fitting
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod0 = a0.integrate()

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################

###### fig set up
fig3, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6)) #fig creationg of 1 by 2
fig3.suptitle('Pro in 0 H Model') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.85, wspace = 0.25, hspace = 0.30)

ax0.semilogy()
ax0.set_title('Pro dynamics ')
ax1.set_title('Model comparison to data')

ax0.set_xlabel('days')
ax0.set_ylabel('cell concentration')
ax1.set_ylabel('Final P value')
ax1.set_xlabel('iteration number ')

l3 = ax0.legend(loc = 'lower right')
l3.draw_frame(False)


#graphing data from df to see 2 different biological reps represented

ax0.errorbar(df0['time'],df0['abundance'],yerr=df0['std1'], marker='o', label = 'avg1')
ax0.errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='o', label = 'avg2')

ax0.plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' model best fit')
a0.plot_uncertainty(ax0,posteriors0,'P',100)
#ax1.scatter(a0res['res'], a0res['abundance'],label = '0H case')
#printing off graph
plt.show()




#########################################################
#graphing model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,5,figsize=[20,7])
#set titles and config graph 
fig4.suptitle('Monoculture parameters in 0 HOOH ')
ax4[0].set_title('Model Dynamic output')
ax4[1].set_title('P0')
ax4[2].set_title('N0')
ax4[3].set_title('k1')
ax4[4].set_title('k2')

ax4[2].set_xlabel('Parameter Value Frequency', fontsize = 16)
#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
#shift fig subplots
fig4.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)


#graph data, model, and uncertainty 
ax4[0].plot(df0.time, df0.abundance, marker='o',label = 'Pro Mono - 0 H ')
ax4[0].plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' Model P best fit')
a0.plot_uncertainty(ax4[0],posteriors0,'P',100)

ax0.plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' Model P best fit')
a0.plot_uncertainty(ax0,posteriors0,'P',100)


# plot histograms of parameter search results 
ax4[1].hist(posteriors0.P0)
ax4[2].hist(posteriors0.N0)
ax4[3].hist(posteriors0.k1)
ax4[4].hist(posteriors0.k2)

#show full graph 
plt.show()

pframe = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe.to_csv('../data/inits/pro9215_inits0.csv')

fig4.savefig('../figures/pro_odelib0')

print("I'm done bro! ")




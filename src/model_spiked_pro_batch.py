
'''

name:   model_spiked_pro_batch.py 

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
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

df = df_mono
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['log4'] = np.log(df['rep4'])
df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
df['avg2'] = df[['rep2', 'rep4']].mean(axis=1)
df['abundance'] = df[['rep1','rep2','rep3', 'rep4']].mean(axis=1)
df['std1'] = df[['rep1', 'rep3']].std(axis=1)
df['std2'] = df[['rep2', 'rep4']].std(axis=1)
df['sigma'] = df[['rep1','rep2','rep3', 'rep4']].std(axis=1)

df['lavg1'] = df[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
df['lavg2'] = df[['log2', 'log4']].mean(axis=1)
df['log_abundance'] = df[['log1','log2', 'log3','log4']].mean(axis=1)
df['stdlog1'] = df[['log1', 'log3']].std(axis=1) #taking stdv of logged reps
df['stdlog2'] = df[['log2', 'log4']].std(axis=1)
df['log_sigma'] = df[['log1','log2', 'log3','log4']].std(axis=1)



#setting working df for model as far as abundance and log abundance values 
#df.rename(columns = {'avg1':'abundance'}, inplace = True) #reaneme main df column to be fit by odelib 
#df.rename(columns = {'lavg1':'log_abundance'}, inplace = True) #reaneme log of main df column to be fit by odelib 

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False) & (df['Vol_number']== 1)]  #assay 0 H 
df4 = df.loc[(df['assay'].str.contains('4', case=False)) & (df['Vol_number']== 1)]


df = df4

#####################################################
#config data in df from raw for odelib usefulness
#####################################################

#making avg columns of technical reps (std here only for graphing, not logged here)

#splitting df of Pro into 0 and 400 H assays 
 

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig2.suptitle('Pro  Monoculture in 400 nM HOOH')
fig2.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)


#format fig  
ax0.set_title('Pro dynamics') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('HOOH dynamics ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('HOOH (nM)')
#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['avg1'],yerr=df4[df4['organism']=='P']['std1'], marker='o', label = 'avg1')
ax0.errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['avg2'],yerr=df4[df4['organism']=='P']['std2'], marker='v', label = 'avg2 ')
ax0.errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['abundance'],yerr=df4[df4['organism']=='P']['sigma'], marker='d', label = 'MEAN ')
# graph 400 H assay even and odd avgs
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['avg1'],yerr=df4[df4['organism']=='H']['std1'], marker='o', label = 'avg1')
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['avg2'],yerr=df4[df4['organism']=='H']['std2'], marker='v', label = 'avg2')
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['sigma'], marker='d', label = 'MEAN')

plt.legend()

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits4 = pd.read_csv("../data/inits/pro9215_inits4.csv")

#setting how many MCMC chains you will run 
nits = 10000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['P','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.000002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.6})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
#setting state variiable  prior guess
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+6})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':100})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
P0_mean = 50000
N0_mean = 900000
H0_mean = 1000

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','P0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          P0 = P0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          P = P0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi = params[0], params[1], params[2],params[3]
    P,N,H = y[0],y[1],y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P - kdam*P*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    dHdt = - phi*P*H
    return [dPdt,dNdt,dHdt]

#df0.loc[:,'log_abundance'] = np.log(10**df0.log_abundance)

# get_models
a4 = get_model(df4) 

#broken here!!!!!!!!!!
# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,static_parameters =set(['k1','k2','N0']),print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################

###### fig set up
fig3, ax3 = plt.subplots(1,2,figsize = (7,6)) #fig creationg of 1 by 2
fig3.suptitle('Pro in 400 H Model') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.85, wspace = 0.25, hspace = 0.30)

ax3[0].semilogy()
ax3[0].set_title('Pro dynamics ')
ax3[1].set_title('HOOH dynamics')

ax3[0].set_xlabel('days')
ax3[0].set_ylabel('cell concentration')
ax3[1].set_ylabel('HOOH concentration')
ax3[1].set_xlabel('days ')

l3 = ax3[0].legend(loc = 'lower right')
l3.draw_frame(False)


#graphing data from df to see 2 different biological reps represented

ax3[0].errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['abundance'],yerr=df4[df4['organism']=='P']['std1'], marker='o', label = 'Mean P')
ax3[0].plot(mod4.time,mod4['P'],color ='r',lw=1.5,label=' P model best fit')
a4.plot_uncertainty(ax3[0],posteriors4,'P',100)

ax3[1].errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['std1'], marker='o', label = 'Mean H')
ax3[1].plot(mod4.time,mod4['H'],color ='r',lw=2.0,label=' H model best fit')
a4.plot_uncertainty(ax3[1],posteriors4,'H',100)

#ax1.scatter(a0res['res'], a0res['abundance'],label = '0H case')
#printing off graph
plt.show()




#########################################################
#graphing P model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,4,figsize=[12,8])
#set titles and config graph 
fig4.suptitle('Pro Monoculture parameters in 400 HOOH ')
ax4[0].set_title('Pro Model Dynamic output')
ax4[1].set_title('P0')
ax4[2].set_title('kdam')
ax4[3].set_title('phi')

ax4[2].set_xlabel('Parameter Value Frequency', fontsize = 16)
#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
#shift fig subplots
fig4.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)


#graph data, model, and uncertainty 
ax4[0].plot(df4[df4['organism']=='P']['time'], df4[df4['organism']=='P']['abundance'], marker='o',label = 'Pro Mono - 4 H ')
ax4[0].plot(mod4.time,mod4['P'],color='r',lw=1.5,label=' Model P best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'P',100)

#ax0.plot(mod4.time,mod4['P'],color='r',lw=1.5,label=' Model P best fit')
#a4.plot_uncertainty(a4,posteriors4,'P',10)



# plot histograms of parameter search results 
ax4[1].hist(posteriors4.P0)
ax4[2].hist(posteriors4.kdam)
ax4[3].hist(posteriors4.phi)

#show full graph 
plt.show()
fig4.savefig('../figures/pro_odelib4_params')


#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 
fig5,ax5 = plt.subplots(1,4,figsize=[12,8])
fig5.suptitle('HOOH parmaters ')
ax5[0].set_title('HOOH Model Dynamic output')
ax5[1].set_title('H0')
ax5[2].set_title('kdam')
ax5[3].set_title('phi')

ax5[2].set_xlabel('Parameter Value Frequency', fontsize = 16)
#make legends
l5 = ax5[0].legend(loc = 'upper left')
l5.draw_frame(False)

ax5[0].plot(df4[df4['organism']=='H']['time'], df4[df4['organism']=='H']['abundance'], marker='o',label = 'H data ')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label=' Model P best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)


# plot histograms of parameter search results 
ax5[1].hist(posteriors4.H0)
ax5[2].hist(posteriors4.phi)
ax5[3].hist(posteriors4.kdam)


#show full graph 
plt.show()



pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/pro9215_inits4.csv')


print("I'm done bro! ")




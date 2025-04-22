
'''

name:   model_detoxers_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/HOOH_dynamics/src'
    
author: DKM

goal: Loop model and graph 0 H Pro assay and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop for 0 and 400 using different init files? (need H connected to 0 H first) 

'''

#read in needed packages 
import helpers as hp
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

# get data and visualize uncertainty

df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_2-5-dataset', header = 1)

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'


dfw = df_all

dfw['log1'] = np.log(dfw['rep1'])
dfw['log2'] = np.log(dfw['rep2'])
dfw['log3'] = np.log(dfw['rep3'])
dfw['log4'] = np.log(dfw['rep4'])
dfw['avg1'] = dfw[['rep1', 'rep3']].mean(axis=1)
dfw['avg2'] = dfw[['rep2', 'rep4']].mean(axis=1)
dfw['abundance'] = dfw[['rep1','rep2','rep3', 'rep4']].mean(axis=1)
dfw['std1'] = dfw[['rep1', 'rep3']].std(axis=1)
dfw['std2'] = dfw[['rep2', 'rep4']].std(axis=1)
dfw['sigma'] = dfw[['rep1','rep2','rep3', 'rep4']].std(axis=1)

dfw['lavg1'] = dfw[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
dfw['lavg2'] = dfw[['log2', 'log4']].mean(axis=1)
dfw['log_abundance'] = dfw[['log1','log2', 'log3','log4']].mean(axis=1)
dfw['stdlog1'] = dfw[['log1', 'log3']].std(axis=1) #taking stdv of logged reps
dfw['stdlog2'] = dfw[['log2', 'log4']].std(axis=1)
dfw['log_sigma'] = dfw[['log1','log2', 'log3','log4']].std(axis=1)

dfw['log_sigma'] = 0.2
dfw.loc[dfw['organism'] == 'H', 'log_sigma'] = 0.08
#######
df_P = df_all[df_all['organism']== 'P']
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  

#####################
#strain slice through vol number selection 
#57,58,59,60
##########################################

vol = 58
#vol58 c
c0 = 'blueviolet'
c1 = 'pink'

df = df_mono.loc[(df_mono['Vol_number'] == vol )].copy()  
df = df.loc[(df['organism'] == 'D' )].copy()  

vol = int(sys.argv[1])
if vol == 57:
    #vol57 colors S
    c0 = 'violet'
    c1 = 'crimson'
    mylab = '$Micromonas$ $commoda$'
elif vol == 59:
    #vol59 colors
    c0 = 'mediumorchid'
    c1 = 'lightcoral'
    mylab = '$Micromonas$ $pusilla$'
elif vol == 58:
    #vol58 c
    c0 = 'blueviolet'
    c1 = 'pink'
    mylab = '$Ostreococcus$ $lucimarinus$'
elif vol == 60:
    #vol60 colors 
    c0 = 'mediumpurple'
    c1 = 'magenta'
    mylab = '$Ostreococcus$ $tauri$'

df0 = df.loc[~df['assay'].str.contains('4', case=False)].copy()
#####################################################
#config data in df from raw for odelib usefulness
#####################################################

#making avg columns of technical reps (std here only for graphing, not logged here)
#splitting df of Pro into 0 and 400 H assays 

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0)= plt.subplots(1,1,figsize = (10,6))
fig2.suptitle('Het Vol '+str(vol)+'  Monoculture in 0 nM HOOH')
fig2.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)

#format fig  
ax0.set_title('Het Vol '+str(vol)+' dynamics') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data

ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 

#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df0[df0['organism']=='D']['time'],df0[df0['organism']=='D']['avg1'],yerr=df0[df0['organism']=='D']['std1'],color = 'b', marker='o', label = 'avg1')
ax0.errorbar(df0[df0['organism']=='D']['time'],df0[df0['organism']=='D']['avg2'],yerr=df0[df0['organism']=='D']['std2'], color = 'g',marker='v', label = 'avg2 ')
ax0.errorbar(df0[df0['organism']=='D']['time'],df0[df0['organism']=='D']['abundance'],yerr=df0[df0['organism']=='D']['sigma'], color = c0,marker='d', label = 'MEAN ')

plt.legend()

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits0 = pd.read_csv('../data/inits/Het_'+str(vol)+ '_inits0.csv')

#setting how many MCMC chains you will run 
nits = 10000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['D','N'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search

#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.00002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.6})
D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})

#still not sure what part of fitting algor this is used for
D0_mean = inits0['D0'][0]
N0_mean = inits0['N0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','D0','N0',],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          D0 = D0_prior.copy(),
                          N0  = N0_prior.copy(),
                          t_steps=1000,
                          D = D0_mean,
                          N = N0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2 = params[0], params[1]
    D,N = max(y[0],0),max(y[1],0)
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dDdt = (k2 * N /( (ksp) + N) )*D 
    dNdt =  - (k2 * N /( (ksp) + N) )*D
    return [dDdt,dNdt]

def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)


# get_models
a4 = get_model(df0) 


# do fitting
posteriors4 = a4.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod4 = a4.integrate()

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################


#########################################################
#graphing P model vs data and params histograms 
#########################################################

fig4,ax4 = plt.subplots(1,3,figsize=[12,4])

fig4.subplots_adjust(wspace=0.3,hspace=0.3)

ax4[0].set_xlabel('Time (days)', fontsize = 12)
ax4[1].set_xlabel('Initial cell density \n($P_{i,0}$, Cells mL$^{-1}$)', fontsize = 12)
ax4[1].set_ylabel('Probability density (x10$^{-5}$)', fontsize = 12)
ax4[2].set_xlabel('Growth rate ($\mu_i$, day$^{-1}$)', fontsize = 12)
ax4[2].set_ylabel('Probability density', fontsize = 12)
ax4[0].set_xlabel('Time (days)', fontsize = 14)
ax4[0].set_ylabel('Cells (mL$^{-1}$)', fontsize = 12)
ax4[0].tick_params(axis='x', labelsize=12)
ax4[0].tick_params(axis='y', labelsize=12)
ax4[1].tick_params(axis='x', labelsize=12)
ax4[1].tick_params(axis='y', labelsize=12)
ax4[2].tick_params(axis='x', labelsize=12)
ax4[2].tick_params(axis='y', labelsize=12)

fig4.subplots_adjust(wspace = 0.3, bottom=0.2)

#shift fig subplots

#graph data, model, and uncertainty 
ax4[0].plot(df0[df0['organism']=='D']['time'], df0[df0['organism']=='D']['abundance'], color = c0, marker='o')
ax4[0].errorbar(df0[df0['organism']== 'D']['time'], df0[df0['organism']== 'D']['abundance'], yerr = df0[df0['organism']== 'D']['sigma'], color = c0, marker='o',label = mylab)
ax4[0].plot(mod4.time,mod4['D'],color='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'D',100)

# plot histograms of parameter search results 
ax4[1].hist(posteriors4.D0,density=True, facecolor = c0)
ax4[2].hist(posteriors4.k2,density=True, facecolor = c0)

#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
#show full graph 

# rescale 
ticks = ax4[1].get_yticks()
ax4[1].set_yticklabels([f"{tick * 1e+5:.0f}" for tick in ticks])

pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/Het_'+str(vol)+ '_inits0.csv')

ax4[0].semilogy()

#show full graph
fig4.savefig('../figures/Het-'+str(vol)+ '_no_spike',bbox_inches='tight')



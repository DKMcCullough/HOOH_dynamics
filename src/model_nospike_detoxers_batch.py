
'''

name:   model_spiked_detoxers_batch.py 

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
import statsmodels.formula.api as smf

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

#######
df_P = df_all[df_all['organism']== 'P']
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  

#####################

#strain slice through vol number selection 
#57,58,59,60
##########################################

vol = int(sys.argv[1])

df = df_mono.loc[(df_mono['Vol_number'] == vol )].copy()  

#vol58 c
c0 = 'blueviolet'
c1 = 'pink'

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

df4 = df.loc[~df['assay'].str.contains('4', case=False)].copy()
#df4 = df4[df4.time < 3] #to keep growth bump in last days to thro off death (kdma) range.

df4['log_sigma'] = hp.get_uncertainty(df4[df4.organism=='D'])
df4.loc[df4['organism'] == 'H', 'log_sigma'] = hp.get_uncertainty(df4[df4.organism=='H'])

if vol == 60:
    df4 = df4[df4.time <5]

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits4 = pd.read_csv('../data/inits/Het_'+str(vol)+ 'nospike_inits4.csv')

# use OLS to constrain growth rate
dfregress = df4[df4.organism=='D']
model = smf.ols('log_abundance ~ time', data=dfregress).fit()
inits4['k2'] = model.params.time

#setting how many MCMC chains you will run 
nits = 100000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['D','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search

#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.00002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':inits4['k2'][0]})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.006})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':12})
#setting state variiable  prior guess
D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':350})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
D0_mean = inits4['D0'][0]
N0_mean = inits4['N0'][0]
H0_mean = inits4['H0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','D0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
                          D0 = D0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          D = D0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    D,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dDdt = (k2 * N /( (ksp) + N) )*D - kdam*D*H
    dNdt =  - (k2 * N /( (ksp) + N) )*D
    dHdt = 12- phi*D*H
    return [dDdt,dNdt,dHdt]

# get_models
a4 = get_model(df4) 

# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,static_parameters =set(['k1','k2','N0']),print_report=True) #, )

# run model with optimal params
mod4 = a4.integrate()

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################


#########################################################
#graphing P model vs data and params histograms 
#########################################################

figall,axall = plt.subplots(2,3,figsize=[12,8])
figall.subplots_adjust(wspace=0.3,hspace=0.3)
ax4,ax5 = axall[0,:],axall[1,:]
#figall.suptitle('Het Vol '+str(vol)+'  Monoculture in 400 nM HOOH')
# set up graph
#set titles and config graph 

ax4[0].set_ylabel('Cells (mL$^{-1}$)', fontsize = 12)
ax4[0].set_xlabel('Time (days)', fontsize = 12)
ax4[1].set_xlabel('log$_{10}$ Initial cell density \n ($P_{i,0}$, cells mL$^{-1}$)', fontsize = 12)
ax4[1].set_ylabel('Probability density', fontsize = 12)
ax4[2].set_xlabel('log$_{10}$ Damage rate \n ($\kappa_{dam,i}$, mL pmol$^{-1}$ day$^{-1}$)', fontsize = 12)
ax4[2].set_ylabel('Probability density', fontsize = 12)

ax4[1].tick_params(axis='x', labelsize=14)
ax4[1].tick_params(axis='y', labelsize=14)
ax4[2].tick_params(axis='x', labelsize=14)
ax4[2].tick_params(axis='y', labelsize=14)

#shift fig subplots

#graph data, model, and uncertainty 
ax4[0].plot(df4[df4['organism']=='D']['time'], df4[df4['organism']=='D']['abundance'], color = c0, marker='o')
ax4[0].errorbar(df4[df4['organism']== 'D']['time'], df4[df4['organism']== 'D']['abundance'], yerr = df4[df4['organism']== 'D']['sigma'], color = c0, marker='o',label = mylab)
ax4[0].plot(mod4.time,mod4['D'],color='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'D',100)

# plot histograms of parameter search results 
hp.sns.kdeplot(np.log10(posteriors4.D0),color=c0,fill=True,ax=ax4[1],linewidth=3)
hp.sns.kdeplot(np.log10(posteriors4.kdam),color=c0,fill=True,ax=ax4[2],linewidth=3)

#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
#show full graph 

#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 
ax5[0].set_ylabel(r'H$_2$O$_2$ concentration (pmol mL$^{-1}$)', fontsize=12)
ax5[0].set_xlabel('Time (Days)', fontsize = 12)
ax5[1].set_xlabel('log$_{10}$ Initial H$_2$O$_2$ concentration \n ($H_0$, pmol mL$^{-1}$)', fontsize = 12)
ax5[1].set_ylabel('Probability density', fontsize = 12)
ax5[2].set_xlabel('log$_{10}$ Detoxification rate \n ($\phi_{det,i}$, pmol cell$^{-1}$ day$^{-1}$)', fontsize = 12)
ax5[2].set_ylabel('Probability density', fontsize = 12)

ax5[1].tick_params(axis='x', labelsize=14)
ax5[1].tick_params(axis='y', labelsize=14)
ax5[2].tick_params(axis='x', labelsize=14)
ax5[2].tick_params(axis='y', labelsize=14)

#plot dynamics and models
ax5[0].plot(df4[df4['organism']=='H']['time'], df4[df4['organism']=='H']['abundance'], color = c1, marker='o')
ax5[0].errorbar(df4[df4['organism']== 'H']['time'], df4[df4['organism']== 'H']['abundance'], yerr = df4[df4['organism']== 'H']['sigma'], color = c0, marker='o',label = 'Hydrogen peroxide')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label='Model best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)

# plot histograms of parameter search results 
hp.sns.kdeplot(np.log10(posteriors4.H0),color=c1,fill=True,ax=ax5[1],linewidth=3)
hp.sns.kdeplot(np.log10(posteriors4.phi),color=c1,fill=True,ax=ax5[2],linewidth=3)

#make legends
l5 = ax5[0].legend(loc = 'upper right')
l5.draw_frame(False)
#show full graph 

for (ax,l) in zip(axall.flatten(),'abcdef'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

ax4[0].semilogy()
ax5[0].semilogy()

figall.subplots_adjust(wspace=0.4,hspace=0.4)
figall.savefig('../figures/Het-'+str(vol)+'nospike_Hparams',bbox_inches='tight')

pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/Het_'+str(vol)+ 'nospike_inits4.csv')

posteriors4.to_csv('../data/posteriors/pos_het_'+str(vol)+ 'nospike.csv')

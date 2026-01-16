
'''

name:   model_spiked_syn_batch.py 

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
df = hp.get_data('coculture')

vol = 53
c0 = 'steelblue'
c1 = 'chocolate'

vol = int(sys.argv[1])
if vol == 52:
    #vol52 colors Syn WH7803
    c0 = 'cornflowerblue'
    c1 = 'darkorange'
    figlab = '$Synechococcus$ WH7803'
elif vol == 28:
    #vol28 colors Syn WH7803 also 
    c0 = 'dodgerblue'
    c1 = 'tomato'
    figlab = '$Synechococcus$ WH7803'
elif vol == 53:
    #vol53 colors WSyn CC9605 
    c0 = 'steelblue'
    c1 = 'chocolate'
    figlab = '$Synechococcus$ CC9605'
elif vol == 54:
    #vol54 colors Syn WH7802 
    c0 = 'darkcyan'
    c1 = 'lightcoral'
    figlab = '$Synechococcus$ WH7802'

#slicing data into abiotic, biotic, and Pro only dataframes
df4 = df.loc[(df['assay'].str.contains('4', case=False)) & (df['Vol_number']== vol)]
sigma4H = hp.get_uncertainty(df4[df4.organism=='H'])
sigma4S = hp.get_uncertainty(df4[df4.organism=='S'])
df4.loc[df['organism'] == 'H', 'log_sigma'] = sigma4H
df4.loc[df['organism'] == 'S', 'log_sigma'] = sigma4S
figa,axa = hp.plot_uncertainty(df4[df4.organism=='S'],sigma4S)
#figa.savefig('../figures/error_syn')

df = df4

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits4 = pd.read_csv('../data/inits/syn_vol'+str(vol)+ '_inits4.csv')
inits0 = pd.read_csv('../data/inits/syn_vol'+str(vol)+ '_inits0.csv')

#setting how many MCMC chains you will run 
nits = 100000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['S','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search

#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.00002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':inits0['k2'][0]})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':12})
#setting state variiable  prior guess
S0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':100})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
S0_mean = inits4['S0'][0]
N0_mean = inits4['N0'][0]
H0_mean = inits4['H0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','S0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
                          S0 = S0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          S = S0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    S,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dSdt = (k2 * N /( (ksp) + N) )*S - kdam*S*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*S
    dHdt = 12- phi*S*H
    return [dSdt,dNdt,dHdt]

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

#########################################################
#graphing P model vs data and params histograms 
#########################################################

figall,axall = plt.subplots(2,3,figsize=[12,8])
figall.subplots_adjust(wspace=0.3,hspace=0.3)
ax4,ax5 = axall[0,:],axall[1,:]

# set up graph
#set titles and config graph 
ax4[0].semilogy()

ax4[0].set_ylabel('Cells (mL$^{-1}$)', fontsize = 12)
ax4[0].set_xlabel('Time (days)', fontsize = 12)
ax4[1].set_xlabel('log$_{10}$ Initial cell density \n ($P_{i,0}$, cells mL$^{-1}$)', fontsize = 12)
ax4[1].set_ylabel('Probability density', fontsize = 12)
ax4[2].set_xlabel('log$_{10}$ Damage rate \n ($\kappa_{dam,i}$, mL pmol$^{-1}$ day$^{-1}$)', fontsize = 12)
ax4[2].set_ylabel('Probability density', fontsize = 12)

ax4[1].tick_params(axis='x', labelsize=12)
ax4[1].tick_params(axis='y', labelsize=12)
ax4[2].tick_params(axis='x', labelsize=12)
ax4[2].tick_params(axis='y', labelsize=12)

#shift fig subplots

#graph data, model, and uncertainty 
ax4[0].plot(df4[df4['organism']=='S']['time'], df4[df4['organism']=='S']['abundance'], color = c0, marker='o',label = figlab)
ax4[0].plot(mod4.time,mod4['S'],color='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'S',100)

# plot histograms of parameter search results 
hp.sns.kdeplot(np.log10(posteriors4.S0),color=c0,fill=True,ax=ax4[1],linewidth=3)
hp.sns.kdeplot(np.log10(posteriors4.kdam),color=c0,fill=True,ax=ax4[2],linewidth=3)

#make legends
l4 = ax4[0].legend(loc = 'lower right')
l4.draw_frame(False)
#show full graph 

#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 
ax5[0].semilogy()

#HOOH dynamics 
ax5[0].set_ylabel(r'H$_2$O$_2$ concentration (pmol mL$^{-1}$)', fontsize=12)
ax5[0].set_xlabel('Time (Days)', fontsize = 12)
ax5[1].set_xlabel('log$_{10}$ Initial H$_2$O$_2$ concentration \n ($H_0$, pmol mL$^{-1}$)', fontsize = 12)
ax5[1].set_ylabel('Probability density', fontsize = 12)
ax5[2].set_xlabel('log$_{10}$ Detoxification rate \n ($\phi_{det,i}$, pmol cell$^{-1}$ day$^{-1}$)', fontsize = 12)
ax5[2].set_ylabel('Probability density', fontsize = 12)

ax5[1].tick_params(axis='x', labelsize=12)
ax5[1].tick_params(axis='y', labelsize=12)
ax5[2].tick_params(axis='x', labelsize=12)
ax5[2].tick_params(axis='y', labelsize=12)

#plot dynamics and models
ax5[0].plot(df4[df4['organism']=='H']['time'], df4[df4['organism']=='H']['abundance'], color = c1, marker='o',label = 'Hydrogen peroxide')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)

# plot histograms of parameter search results 
hp.sns.kdeplot(np.log10(posteriors4.H0),color=c1,fill=True,ax=ax5[1],linewidth=3)
hp.sns.kdeplot(np.log10(posteriors4.phi),color=c1,fill=True,ax=ax5[2],linewidth=3)

#make legends
l5 = ax5[0].legend(loc = 'lower left')
l5.draw_frame(False)
#show full graph 

for (ax,l) in zip(axall.flatten(),'abcdef'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

figall.subplots_adjust(wspace=0.35,hspace=0.35)
figall.savefig('../figures/syn_'+str(vol)+'_Hparams',bbox_inches='tight')
figall.savefig('../figures/syn_'+str(vol)+'_Hparams.tiff',bbox_inches='tight',format='tiff',dpi=200)

pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/syn_vol'+str(vol)+ '_inits4.csv')

posteriors4.to_csv('../data/posteriors/pos_syn_'+str(vol)+ 'spike.csv')

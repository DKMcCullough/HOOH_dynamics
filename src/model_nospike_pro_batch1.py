
'''

name:   model_spiked_pro_batch1.py 

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

#slicing data into abiotic, biotic, and Pro only dataframes
df4 = df.loc[(~df['assay'].str.contains('4', case=False)) & (df['Vol_number']== 1)]
df4 = df4[df4.time < 6] #to keep growth bump in last days to thro off death (kdma) range.

df = df4

c0 = 'g'
c1 = 'darkgoldenrod'

#slicing data into abiotic, biotic, and Pro only dataframes
#sigma4H = hp.get_uncertainty(df4[df4.organism=='H'])
#sigma4P = hp.get_uncertainty(df4[df4.organism=='P'])
#df4.loc[df['organism'] == 'H', 'log_sigma'] = sigma4H
#df4.loc[df['organism'] == 'P', 'log_sigma'] = sigma4P
#figa,axa = hp.plot_uncertainty(df4[df4.organism=='P'],sigma4P)
#figa.savefig('../figures/error_pro')

df = df4


#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits0 = pd.read_csv("../data/inits/pro_MIT9215_inits0_1.csv")

#setting how many MCMC chains you will run 
nits = 100000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['P','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':2e+5})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':inits0['k2'][0]})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':inits0['kdam'][0]})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':inits0['phi'][0]})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':100})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
P0_mean = inits0['P0'][0]
N0_mean = inits0['N0'][0]
H0_mean = inits0['H0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','P0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
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
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    P,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P - kdam*P*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    dHdt = 10 - phi*P*H
    return [dPdt,dNdt,dHdt]


#####################################################
#modeling and fitting 
#####################################################

# get_models
a4 = get_model(df4) 

# do fitting
posteriors4 = a4.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,static_parameters =set(['k1','k2','N0']),print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()

#########################################################
#graphing P model vs data and params histograms 
#########################################################

# set up graph
figall,axall = plt.subplots(2,3,figsize=[12,8])
figall.subplots_adjust(wspace=0.3,hspace=0.3)
ax4,ax5 = axall[0,:],axall[1,:]
#set titles and config graph 

ax4[0].set_ylabel('Cells (mL$^{-1}$)', fontsize=12)
ax4[0].set_xlabel('Time (days)',fontsize=12)

ax4[1].set_xlabel('Initial cell density ($P_{i,0}$, cells mL$^{-1}$)', fontsize = 12)
ax4[1].set_ylabel('Probability density (x10$^{-6}$)', fontsize = 12)
ax4[2].set_xlabel('Damage rate \n ($\kappa_{dam,i}$, x10$^{-5}$ mL pmol$^{-1}$ day$^{-1}$)', fontsize = 12)
ax4[2].set_ylabel('Probability denity', fontsize = 12)

#graph data, model, and uncertainty 
ax4[0].plot(df4[df4['organism']=='P']['time'], df4[df4['organism']=='P']['abundance'], color =  c0,marker='o')
ax4[0].errorbar(df4[df4['organism']== 'P']['time'], df4[df4['organism']== 'P']['abundance'], yerr = df4[df4['organism']== 'P']['sigma'], color = c0, marker='o',label = r'$Prochlorococcus$')
a4.plot_uncertainty(ax4[0],posteriors4,'P',100)
ax4[0].plot(mod4.time,mod4['P'],color='r',lw=1.5,label='Model best fit')

# plot histograms of parameter search results 
ax4[1].hist(posteriors4.P0, color =  c0,density=True)
ax4[2].hist(posteriors4.kdam*1e+5,color = c0,density=True)

# rescale 
ticks = ax4[1].get_yticks()
ax4[1].set_yticklabels([f"{tick * 1e+6:.0f}" for tick in ticks])

#make legends
l4 = ax4[0].legend(loc = 'lower left')
l4.draw_frame(False)
#show full graph 

#ax4[0].set_ylim([100, 5000000])
ax4[0].semilogy()

for (ax,l) in zip(axall.flatten(),'abcdef'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#########################################################
#graphing H model vs data and params histograms 
#########################################################

#ax5[0].set_ylim([200, 500])

#HOOH dynamics 
ax5[0].set_ylabel(r'H$_2$O$_2$ concentration (pmol mL$^{-1}$)', fontsize=12)
ax5[0].set_xlabel('Time (days)',fontsize=12)
ax5[1].set_xlabel('Initial H$_2$O$_2$ concentration ($H_0$, pmol mL$^{-1}$)', fontsize = 12)
ax5[1].set_ylabel('Probability density', fontsize = 12)
ax5[2].set_xlabel('Detoxification rate \n ($\phi_{det,i}$, x10$^{-6}$ pmol cell$^{-1}$ day$^{-1}$)', fontsize = 12)
ax5[2].set_ylabel('Probability density', fontsize = 12)

ax5[0].plot(df4[df4['organism']=='H']['time'], df4[df4['organism']=='H']['abundance'], color =  c1,marker='o')
ax5[0].errorbar(df4[df4['organism']== 'H']['time'], df4[df4['organism']== 'H']['abundance'], yerr = df4[df4['organism']== 'H']['sigma'], color = c0, marker='o',label = 'Hydrogen peroxide')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label='Model best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)

l5 = ax5[0].legend(loc = 'lower left')
l5.draw_frame(False)

# plot histograms of parameter search results 
ax5[1].hist(posteriors4.H0,color =  c1, density=True)
ax5[2].hist(posteriors4.phi*1e+6, color = c1,density=True)


figall.savefig('../figures/pro1_0nm_h2o2',bbox_inches='tight')
figall.savefig('../figures/figure5test.tiff',dpi=200,format='tiff',bbox_inches='tight')

pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/pro_MIT9215_inits0_1.csv')

posteriors4.to_csv('../data/posteriors/pos_pro_no_spike.csv')



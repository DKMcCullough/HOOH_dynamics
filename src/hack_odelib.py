#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:35:34 2023
name : hack_odelib

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/Project_6_BCCcocultures/src'
@author: dkm
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy.integrate import odeint
import ODElib
import random as rd

################
#data import
#################


#main_dir = "/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/Project_6_BCCcocultures/"

df_all = pd.read_csv("odelib_hack_dt_emily.csv")
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df = df_all
df['log_sigma'] = df['log_sigma'].clip(lower=0) 
df = df.rename(columns={"raw_abundance": "abundance"})

inits = pd.read_csv("inits/abiotic.csv")

df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)] 
df_P = df.loc[df['assay'].str.contains('P_', case=False)] 
df_S = df.loc[df['assay'].str.contains('S_', case=False)] 
df_co = df.loc[df['assay'].str.contains('co', case=False)] 

#df_abiotic = df.loc[df['assay'].str.contains("abiotic", case=False)].copy()
#df_Syn = df.loc[~df['assay'].str.contains("abiotic", case=False)].copy()

#S0 = Syn_0.loc[(Syn_0['time'] == 0) & (Syn_0['organism'] == 'S'), 'abundance'].iloc[0]
#H0 = Syn_0.loc[(Syn_0['time'] == 0) & (Syn_0['organism'] == 'H'), 'abundance'].iloc[0]


df0 = df_abiotic.loc[~ df_abiotic['assay'].str.contains('4', case=False)] 
df4 = df_abiotic.loc[df_abiotic['assay'].str.contains('4', case=False)] 




##############
#param set up 
##############


deltah = 0.02       #decay rate of HOOH via Syn 
Sh = 6


##################
#functions
###################



def set_best_params(model,posteriors,snames):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    bestchain = posteriors.iloc[im]["chain#"]
    posteriors = posteriors[posteriors["chain#"]==bestchain]
    model.set_parameters(**posteriors.loc[im][a1.get_pnames()].to_dict())
    model.set_inits(**{o:posteriors.loc[im][a1.get_pnames()].to_dict()[o+'0'] for o in ['H']})

def plot_uncertainty(ax,model,posteriors,ntimes):
    for a in range(ntimes):
        im = rd.choice(posteriors.index)
        model.set_inits(**{'H':posteriors.loc[im][model.get_pnames()].to_dict()['H0']})
        model.set_parameters(**posteriors.loc[im][model.get_pnames()].to_dict())
        mod = model.integrate()
        ax.plot(mod.time,mod['H'],c=str(0.8),lw=1,zorder=1)

# initialize the class for fitting
def get_model(df):
    a1=ODElib.ModelFramework(ODE=abiotic,
                          parameter_names=['deltah','Sh', 'H0'],
                          state_names = snames,
                          dataframe=df,
                          deltah = deltah_prior.copy(),
                          Sh = Sh_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000000,
                          H = H0_mean,
                         )
    return a1


def abiotic(y,t,params):
    deltah,Sh = params[0], params[1]
    H = y[0]
    dHdt = Sh - deltah*H 
    return [dHdt]


pw = 1

deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})

Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':6})

H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})


H0_mean = df0.loc[df0['time'] == 0, 'abundance'].iloc[0]



# state variable names
snames = ['H']

# get_models
a1 = get_model(df0) 
#a4 = get_model(df4) 
 

# define initial conditions
chain_inits = inits

# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 100000

# do fitting
posteriors1 = a1.MCMC(chain_inits=chain_inits,iterations_per_chain=nits,cpu_cores=4)
#posteriors4 = a4.MCMC(chain_inits=chain_inits,iterations_per_chain=nits,cpu_cores=1)

# set best params
set_best_params(a1,posteriors1,snames)
#set_best_params(a4,posteriors4,snames)
# run model with optimal params
mod1 = a1.integrate()
#mod4 = a4.integrate()
#########################################################
# modeling
#########################################################

# pro model graph
f4,ax4 = plt.subplots(1,2,figsize=[8,4.5])
ax4[0].plot(df0.time,df0.abundance, marker='o',label = 'abiotic - 0 H ')
ax4[0].plot(mod1.time,mod1['H'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax4[0],a1,posteriors1,100)

# plot histograms
ax4[1].hist(posteriors1.Sh)
#ax4[1].hist(posteriors4.Sh)


ax4[0].set_xlabel ('Time (days )')
ax4[1].set_xlabel (r'H supply rate, (day$^{-1}$)')
ax4[0].set_ylabel ('Abundance/mL')
ax4[1].set_ylabel ('Frequency')
ax4[0].set_yscale('log')
ax4[0].legend()
plt.show()
'''
# pro model graph
f5,ax5 = plt.subplots(1,2,figsize=[8,4.5])
ax5[0].plot(df4.time,df4.abundance, marker='o',label = 'abiotic - 400 H')
ax5[0].plot(mod1.time,mod1['H'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax4[0],a1,posteriors1,100)

# plot histograms
ax5[1].hist(posteriors1.deltah)
#ax5[1].hist(posteriors4.deltah)


ax5[0].set_xlabel ('Time (days )')
ax5[1].set_xlabel (r'H decay rate, (day$^{-1}$)')
ax5[0].set_ylabel ('Abundance/mL')
ax5[1].set_ylabel ('Frequency')
ax5[0].set_yscale('log')
ax5[0].legend()







a1=ODElib.ModelFramework(ODE=abiotic,
                          parameter_names=['deltah','Sh', 'H0'],
                          state_names = ['H'],
                          dataframe=df0,
                          deltah = deltah_prior.copy(),
                          Sh = Sh_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000000,
                          H = H0_mean,
                         )

chain_inits = inits
#f = a1.plot()
a1posterior = a1.MCMC(chain_inits=chain_inits,iterations_per_chain=500,cpu_cores=1)




a1.plot()
posteriors = a1.MCMC(chain_inits=chain_inits,iterations_per_chain=1000,cpu_cores=1)
im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
bestchain = posteriors.iloc[im]["chain#"]
posteriors = posteriors[posteriors["chain#"]==bestchain]
a1.set_parameters(**posteriors.loc[im][a1.get_pnames()].to_dict())
a1.set_inits(**{o:posteriors.loc[im][a1.get_pnames()].to_dict()[o+'0'] for o in ['H']})
a1.plot()

plt.show()

'''


plt.show()
print('done on the cluster!!!')




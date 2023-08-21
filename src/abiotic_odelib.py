#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:50:19 2022

name:Pro_woHOOH_hack.py

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/Project_6_cocultures/src

author: DKM


goal: ifir base delta and mu for Pro in 0 HOOH 



"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy.integrate import odeint
import ODElib

#SynWH7803 (Vol 28 and 52) is KatG possitive
#SynWH8102 (Vol 54) is KatG negative 


#data import


df_0 = pd.read_csv("Abiotic_0H_ODElib.csv")
df_400 = pd.read_csv('Abiotic_400H_ODElib.csv')

df = df_400
#df_all = df_all.rename({'Time(days)':'times'}, axis=1) 
df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)


#df[df < 0] = 0
df['log_sigma'] = df['log_sigma'].clip(lower=0) #stdv of HOOHs were neg
df.dropna(inplace=True)


#
##################################################3
# parameter and variable Set UP 
#############################################


step = 0.0001
ndays = 8
mtimes = np.linspace(0,ndays,int(ndays/step))


deltah = 0.2       #decay rate of HOOH via Syn 
Sh = 6


H = np.array([])


#initial values to be used for odeint start 

H0 = df['abundance'][0]    #400 actual.....need to get data h0    #nM

#inits = (P0,D0,N0,H0)

#H0 is 400 or 0, SH is 0 both times???
params = (deltah,Sh)
y = [H]
inits = (H0)

def abiotic(y,t,params):
    deltah,Sh = params[0], params[1]
    H = y[0]
    dHdt = Sh - deltah*H 
    return [dHdt]

#lowH = odeint(abiotic,inits, mtimes, args = (params,))
#lHs = lowH[:,0]




deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':1,'scale':1e-7})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':1,'scale':1e-3})

# initialize the class for fitting
a1=ODElib.ModelFramework(ODE=abiotic,
                          parameter_names=['deltah','Sh'], 
                          state_names = ['H'],
                          dataframe=df,
                          deltah = deltah_prior.copy(),
                          Sh = Sh_prior.copy(),
                          t_steps=100000,
                          H = H0,
                         )


f = a1.plot()
a1posterior = a1.MCMC(chain_inits=32,iterations_per_chain=1000,cpu_cores=8,fitsurvey_samples=10000,sd_fitdistance=6.0)

im = a1posterior.loc[a1posterior.chi==min(a1posterior.chi)].index[0]
a1.set_parameters(**a1posterior.loc[im][a1.get_pnames()].to_dict())
#a1.set_inits(**{'H':a1posterior.loc[im][a1.get_pnames()].to_dict()['H0']})

# run the model again, now with fitted parameters
mod = a1.integrate()

# plot fitted model

plt.errorbar(df['time'],df['abundance'], yerr =df['log_sigma'], label = 'data')
plt.plot(a1.times,mod['H'],label='fitted',c='g')
plt.semilogy()
# legend
l = plt.legend()
l.draw_frame(False)

plt.show()
#32,1000,8,10000,6.0







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


df_all = pd.read_csv("Vol1_Vol28_ODElib.csv")

#df_all = df_all.rename({'Time(days)':'times'}, axis=1) 
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

df = df_all
#df[df < 0] = 0
df = df['log_sigma'].clip(lower=0)

'''
#df_all = df_all.iloc[0:32, 0:12]   #need a beter way to do this 
df_abiotic = df_all.loc[df_all['assay'].str.contains("abiotic", case=False)].copy()
df_b = df_all.loc[~df_all['assay'].str.contains("abiotic", case=False)].copy() #b for biotic
#tr = '400' #treatment used...one we use to split between 0 and treatment. 
#df_b['treatment'] = 0  #add treatment column with value of 0 for everything 
#df_b['treatment'].loc[df_b['assay'].str.contains(tr, case=False)] = 400 

#df_pro = df_b[df_b['Vol_number']==1.0] #and get rifd of non vol 28 or 54 cocultures before graphing
#df_pro = df_pro[((~df_pro['coculture']==True))]

#df_pro['abundance'] = np.nanmean(np.r_[[df_pro[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0)
#df_pro['biomass_stdv'] = np.std(np.r_[[df_pro[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0)


#treats = (0,400)

fig2, (ax1,ax2)= plt.subplots(2,1,sharex = True, figsize = (6,7))
fig2.suptitle('Prochlorococcus Monoculture Dynamics')
ax1.set_ylabel('Pro Biomass (cells/mL)')
ax2.set_ylabel('HOOH (nM)')
ax2.set_xlabel('Time (days)')

for t in treats:
    print(t)
    df = df_pro.loc[df_pro['treatment']==int(t)].copy()
    ax1.plot(df['times'],df['abundance'], marker='o' ,label='Pro growth '+ str(t)) #,yerr = df['biomass_stdv'],
    ax2.plot(df['times'], df['HOOH_avg'], marker = 'o',label=str(t) + ' with Pro ')#,color = colors[ai])#, label=h_lab) #,yerr = df['HOOH_stdv']
    print(df)


plt.legend()
plt.show()




#df = df_pro.loc[df_pro['treatment']==0].copy()
#df_4 = df_pro.loc[df_pro['treatment']==400].copy()

fig2, (ax1,ax2)= plt.subplots(2,1,sharex = True, figsize = (6,7))
fig2.suptitle('Prochlorococcus in HOOH ')
ax1.set_ylabel('Pro Biomass (cells/mL)')
ax2.set_ylabel('HOOH (nM)')
ax2.set_xlabel('Time (days)')
ax1.plot(df['times'],df['abundance'], marker='o' ,label='Pro growth no HOOH') #,yerr = df['biomass_stdv'],
ax2.plot(df['times'], df['HOOH_avg'], marker = 'o',label='Pro only ')

'''
#plt.show()

#
##################################################3
# parameter and variable Set UP 
#############################################


step = 0.001
ndays = 8
mtimes = np.linspace(0,ndays,int(ndays/step))


Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 
Qnd = 1#

k1p =  0.00002     #Pro alpha
k2 = 0.83    #************
ksp = k2/k1p
dp = 0.2   #pro delta
kdam = 4.0e-3   #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phip = 1.1e-7    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002
SN = 0
Sh = 8
k1d =  0.000011
dd =  0.2
kdamd = 0.10
ksd = k2/k1d


#empty arrays to be populated by odeint when calling leak function
P = np.array([])
N = np.array([])
H = np.array([])
D  = np.array([])
#y = [P,D,N,H]

#initial values to be used for odeint start 
P0 = 3e5
N0 = 1.0e8        #nM 
H0 = 10    #400 actual.....need to get data h0    #nM
D0 = 5e4
#inits = (P0,D0,N0,H0)

#H0 is 400 or 0, SH is 0 both times???
params = (ksp,k2,dp,deltah,phip,SN,Sh,kdam)
y = [P,N,H]
inits = (P0, N0, H0)

def pro_mono(y,t,params):
    ksp,k2,dp,deltah,phip,SN,Sh,kdam = params[0], params[1], params[2], params[3],params[4], params[5], params[6], params[7]
    P,N,H = y[0],y[1],y[2]
    dPdt = (k2 * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P     
    dNdt =  SN - ((k2 * N /( (ksp) + N) )*P* Qnp) - rho*N    
    dHdt = Sh - deltah*H - phip*H*P #phi being S cell-specific detox rate
    return [dPdt,dNdt,dHdt]

Pro = odeint(pro_mono,inits, mtimes, args = (params,))
Ps = Pro[:,0]
Ns = Pro[:,1]
Hs = Pro[:,2]

#Model graphing 

'''
ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'model P ') #'k1 =' + str(k1p))
ax1.text(3,3e5,s= ('kdam = '+str(kdam)+ ' and phi = '+str(phip)))
ax2.plot(mtimes,Hs,linewidth = 2, color = 'b', label = 'model H')

ax1.semilogy()
ax2.semilogy()
ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
plt.show()



print('*** Done ***')
'''

#odelib attempt 

#df1 = df.loc[:,['times','Plank_type','HOOH_avg','HOOH_stdv','abundance','biomass_stdv']].copy()
#df1['log_sigma'] = np.log(df1['biomass_stdv'])

#organism	time	abundance	log_sigma
#df1 = df1.rename(columns={'times': 'time'})
#df1['organism'] = 0
#df1['organism'].loc[df1['Plank_type'].str.contains('Pro', case=False)] = 'P' 

ksp_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':2,'scale':1e-8})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':1,'scale':0.2})
dp_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':1,'scale':0.2})
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':3,'scale':1e-8})
phip_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':3,'scale':1e-8})
SN_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':1,'scale':1e-3})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':1,'scale':1e-3})
kdam_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':1,'scale':1e-8})

# initialize the class for fitting
P1=ODElib.ModelFramework(ODE=pro_mono,
                          parameter_names=['ksp','k2','dp','deltah','phip','SN','Sh','kdam'], 
                          state_names = ['P','N','H'],
                          dataframe=df,
                          ksp= ksp_prior.copy(),
                          k2 = k2_prior.copy(),
                          dp = dp_prior.copy(),
                          deltah = deltah_prior.copy(),
                          phip = phip_prior.copy(),
                          SN = SN_prior.copy(),
                          Sh = Sh_prior.copy(),
                          kdam = kdam_prior.copy(),
                          t_steps=288,
                          P = P0,
                          H = H0,
                          N = N0
                         )


f = P1.plot()
posterior = P1.MCMC(chain_inits=32,iterations_per_chain=1000,cpu_cores=8,fitsurvey_samples=10000,sd_fitdistance=6.0)


#32,1000,8,10000,6.0







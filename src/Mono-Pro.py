#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:05:08 2023

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src/Monocultures.py'


@author: dkm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint


df_P= pd.read_csv("../data/Mono400_Pros.csv") 
#df_P['abundance'] = np.nanmean(np.r_[[df_P[i] for i in ['rep1','rep2','rep3','rep4']]], axis=0,skipna=True)
#df_P['biomass_stdv'] = np.std(np.r_[[df_P[i] for i in ['rep1','rep2','rep3','rep4']]], axis=0,skipna=True)
df_P['Pavg'] = df_P[['rep1', 'rep2', 'rep3', 'rep4']].mean(axis=1)
df_P['Pstd'] = df_P[['rep1', 'rep2', 'rep3', 'rep4']].std(axis=1)
print(df_P['Pavg'],df_P['Pstd'])
'''
#pros are mit9313, uh18301, mit9312, med4
#syns are wh7803_28, wh_7803_52, cc9605, wh8102
#heterotrophs are rcc299, ccmp1545, ccmp2972A, caoth95
#abiotic is HOOH alone 

'''
#graph params
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 18})

abiotic = pd.read_csv("../data/Mono-abiotic.csv")

strains = df_P['strain'].unique()
nstrains = strains.shape[0]
colors = ['r','orange','g','b','k']
fig1, (ax)= plt.subplots(nstrains,2,figsize = (10,16))

for (S,si) in zip(strains,range(nstrains)):   
    df = df_P[((df_P['strain']==S))].copy()
    df0 = df[((df['assay']=='plus_0'))].copy()
    df400 = df[((df['assay']=='plus_400'))].copy()
    ax[si,0] = df0.plot(kind='scatter', x='Time (days)', y='Pavg', yerr='Pstd',style="-", label = S, title = '0 HOOH assay', ylabel = 'cells per mL',logy = True)
    ax[si,1] = df400.plot(kind='scatter', x='Time (days)', y='Pavg', yerr='Pstd',style="-", label = S, title = '400 HOOH assay', ylabel = 'cells per mL',logy = True)


'''
# define subplot grid
fig, axs = plt.subplots(nrows=nstrains, ncols=2, figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
fig.suptitle("Pro Dynamics", fontsize=18, y=0.95)

# loop through tickers and axes
for S, ax in zip(strains, axs.ravel()):
    # filter df for ticker and plot on specified axes
    df[df["strain"] == S].plot(ax=ax,kind='scatter', x='Time (days)', y='Pavg', yerr='Pstd',style="-", label = S, title = '0 HOOH assay', ylabel = 'cells per mL',logy = True)


plt.show()
'''


#plt.show()

#fig1.suptitle('Pro Dynamics')
#plt.subplots_adjust(wspace = 0.5, top = 0.9,bottom = 0.1)
'''
for S in strains:
    df = df_P[((df_P['strain']==S))].copy()
    df0 = df[((df['assay']=='plus_0'))].copy()
    df400 = df[((df['assay']=='plus_400'))].copy()
    #ax1 = df0.plot(x = 'Time (days)',y = 'Pavg', label = S)
    ax2 = df0.plot(kind='scatter', x='Time (days)', y='Pavg', yerr='Pstd',style="-", label = S, title = '0 HOOH assay', ylabel = 'cells per mL',logy = True)
    ax3 = df400.plot(kind='scatter', x='Time (days)', y='Pavg', yerr='Pstd',style="-", label = S, title = '400 HOOH assay', ylabel = 'cells per mL',logy = True)

'''



'''

#Graph data for each Pro 

df1 = P9313




fig1,(ax1) = plt.subplots(figsize = (7,10))

ax1.set_ylabel('Pro Abundance (cells/mL)')
ax1.set_xlabel('Time (days)')
ax1.plot(neg_growth['times'],neg_growth['abundance'], color = 'goldenrod', marker='o' ,label='Syn WH8102') #,yerr = df['biomass_stdv'],

#plot growth without HOOH 
ax1.plot(neg_growth['times'],neg_growth['abundance'], color = 'goldenrod', marker='o' ,label='Syn WH8102') #,yerr = df['biomass_stdv'],
#ax1.plot(pos_growth['times'],pos_growth['abundance'], marker='o' ,label='Syn WH7803') #,yerr = df['biomass_stdv'],




#model set up 

step = 0.01
ndays = 2.1
times = np.linspace(0,ndays,int(ndays/step))
Qn = (9.4e-15*(1/(14.0))*1e+9)   #Nitrogen Quota for Pro from Bertillison? 

#S_base =   30.0             #3.0 ~160nM from BCC paper calculations 
     
P = (5e5) # units are cells per mL
#S = (S_base + 0)    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022

#k1= 0.002   #ALPHA  
#k2 = 0.22  #VMAX
#nrdelta = 0.06      #nutrient replete delta  #Kdam for this 
#nddelta = 0.19       #nutrient deplete delta kddam for this
#mumax = k2    #set max/min for this known from lit? 

#parameter values (k1 = alpha, k2 = Vmax, kdam = initial HOOH damage, kddam = HOOH damage when N runs out)
k1s = np.r_[[0.02, 0.02, 0.02, 0.02, 0.02]]*1e+100     
k2s = [0.32, 0.32, 0.32, 0.32, 0.32]
kdams = [0.07,0.2, 0.6, 1.3, 3.1]
kddams = [0.03, 0.5, 2, 2, 1.4]

params = list(zip(k1s,k2s,kdams,kddams))

fig2,ax2 = plt.subplots()
#df = dc['df_50']
#data = df['avg_exp']
#S_base =  (data.max()- data.iloc[0])*Qn 

for t in treatments: 
    count = treatments.index(t)
    k1 = params[count][0]
    k2 = params[count][1]
    kdam = params[count][2]
    kddam = params[count][3]
    ks = (k2/k1)   #set max value for this that is know from lit? (between 0.01 and 0.015 for N metabolism in )
    SsEuler = np.array([])
    PsEuler = np.array([])
    P = 1.3e5
    #print(kdam,kddam)
    S_base =   162.0             #3.0 ~160nM from BCC paper calculations 
    df = dc['df_50']
    data = df['avg_exp']
    S =  (data.max()- data.iloc[0])*Qn    #using QN and 50 treatment as base N for test.     
    # units are cells per mL
    #S = S_base    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022
    for t in times:
            PsEuler = np.append(PsEuler,P)
            SsEuler = np.append(SsEuler,S)
            #if (S>2e-3):
            #    delta = kdam
            #else:
            #    delta = kddam
            delta = kdam
            #print(S)
            #print(P)
            dPdt = k2 * P * S /((ks) + S) - delta*P
            dSdt = -P*(k2*(S)/((ks)+S))*Qn
            if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                    S = S + dSdt*step  #S = 4e-47
            else:
                    S = S + dSdt*step
            P = P + dPdt*step
    ax1.plot(times,(PsEuler), linestyle = 'dashed', color = colors[count]) 
    ax2.plot(times,(SsEuler), linestyle = 'dashed', color = colors[count])


    
#ax2_legend = [zip(treatments,colors)]  #ipz doesn't give a printable object for legend to show
ax2.set(xlabel= 'Time (days)', ylabel='Nitrogen (nM)')
ax2.set_title = ('Nutrient dynamics') #not showing up for some reason? 
#ax2.legend([(list(a2_legend))] ,loc = 'lower left')
ax2.semilogy()

'''
print('Done')

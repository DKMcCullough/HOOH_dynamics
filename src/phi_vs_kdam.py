#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:13:24 2023
Created on Tue Nov 15 23:50:19 2022

DDname:phi_vs_kdam.py

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/HOOH_dynamics/src

author: DKM


goal: icompare phi and kdam values 
"""



#SynWH7803 (Vol 28 and 52) is KatG possitive Syn + 
#SynWH8102 (Vol 53) is KatG negative  Syn - 
#Pro = MIT9215


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["font.family"] = "Times New Roman"


df_Pro_1 =  pd.read_csv('../data/inits/pro_MIT9215_inits4_1.csv')
df_Pro_2 = pd.read_csv('../data/inits/pro_MIT9215_inits4_2.csv')
df_Syn_28 = pd.read_csv('../data/inits/syn_vol28_inits4.csv')
df_Syn_52 = pd.read_csv('../data/inits/syn_vol52_inits4.csv')
df_Syn_53 = pd.read_csv('../data/inits/syn_vol53_inits4.csv')
df_Syn_54 = pd.read_csv('../data/inits/syn_vol54_inits4.csv')
#df_Het_55 = pd.read_csv('../data/inits/Het_55_inits4.csv')
df_Het_57 = pd.read_csv('../data/inits/Het_57_inits4.csv')
df_Het_58 = pd.read_csv('../data/inits/Het_58_inits4.csv')
df_Het_59 = pd.read_csv('../data/inits/Het_59_inits4.csv')
df_Het_60 = pd.read_csv('../data/inits/Het_60_inits4.csv')

df_Het_55_3 = pd.read_csv('../data/inits/Het_55_inits_3-spike.csv')
df_Het_55_8 = pd.read_csv('../data/inits/Het_55_inits_8-spike.csv')

df_co_Syn_28 = pd.read_csv('../data/inits/coculture_1,28_pro_inits4.csv')
df_co_Syn_52 = pd.read_csv('../data/inits/coculture_1,52_pro_inits4.csv')
df_co_Syn_53 = pd.read_csv('../data/inits/coculture_1,53_pro_inits4.csv')
df_co_Syn_54 = pd.read_csv('../data/inits/coculture_1,54_pro_inits4.csv')

df_co_Het_57 = pd.read_csv('../data/inits/Het_coculture_1,57_pro_inits4.csv')
df_co_Het_58 = pd.read_csv('../data/inits/Het_coculture_1,58_pro_inits4.csv')
df_co_Het_59 = pd.read_csv('../data/inits/Het_coculture_1,59_pro_inits4.csv')
df_co_Het_60 = pd.read_csv('../data/inits/Het_coculture_1,60_pro_inits4.csv')



df_mono = pd.DataFrame({'names' : ('Syn eWH7803_Vol52','Syn eWH7803_Vol28','Syn eCC9605_Vol53', 'Pro eMIT9215_Vol1_sample1','Pro eMIT9215_Vol1_sample2','Micromonas commoda - Vol_57','Micromonas pusilla - Vol_58','Ostreococcus lucimarinus - Vol_59','Ostreococcus tauri - Vol_60', 'Alteromonas macleodii  - EZ55'),
                                'phis' : [df_Syn_52.phi,df_Syn_28.phi,df_Syn_53.phi,df_Pro_1.phi,df_Pro_2.phi,df_Het_57.phi, df_Het_58.phi,df_Het_59.phi,df_Het_60.phi,df_Het_55_3.phi],
                                'kdams' : [df_Syn_52.kdam,df_Syn_28.kdam,df_Syn_53.kdam,df_Pro_1.kdam,df_Pro_2.kdam,df_Het_57.kdam, df_Het_58.kdam,df_Het_59.kdam,df_Het_60.kdam ,df_Het_55_3.kdam] }, 
                                columns=['names','phis', 'kdams'])


####################################################
#monoculture Pro, Syn, and Het Kdam and Phis
#####################################################
colors = ['cornflowerblue','dodgerblue','steelblue','g','lightgreen','violet','blueviolet','mediumorchid','mediumpurple','b']
markers = ['o','o','o','D','D','v','v','v','v','d'] #HAVE THIS COORS[PMDWOTJ GENES]

#Not all of the data sets have cocurrent HOOH measurements so Phi is less constrained and therefore can be biassing the conclusions a lot. 
fig3, (ax1)= plt.subplots(figsize = (7,6))
fig3.suptitle('Tradeoffs', fontsize = 18)
ax1.set_xlabel('Detoxification rate ($\phi_{max})$',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$)',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

#ax1.set_xlim([1e-11, 1e-5]) #phi
#ax1.set_ylim([1e-5, 1e-2]) #kdam


for i, r in df_mono.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = colors[i], marker = markers[i],markersize=11, linewidth=0.1, label=r['names'])

l3 = plt.legend(loc = 'best', fontsize = 12)

ax1.semilogy()
ax1.semilogx()

fig3.savefig('../figures/tradeoffs')



############################################
#coexistence Phis and Kdams
##############################################


colors = ['cornflowerblue','dodgerblue','steelblue','violet','blueviolet','mediumorchid','mediumpurple']
markers = ['o','o','o','v','v','v','v'] #HAVE THIS COORS[PMDWOTJ GENES]

df_co_help = pd.DataFrame({'names' : ('Syn eWH7803_Vol52','Syn eWH7803_Vol28','Syn eCC9605_Vol53','Micromonas commoda - Vol_57','Micromonas pusilla - Vol_58','Ostreococcus lucimarinus - Vol_59','Ostreococcus tauri - Vol_60'),
                                'phis' : [df_co_Syn_52.phis,df_co_Syn_28.phis,df_co_Syn_53.phis,df_co_Het_57.phis, df_co_Het_58.phis,df_co_Het_59.phis,df_co_Het_60.phis],
                                'kdams' : [ df_co_Syn_52.kdams ,df_co_Syn_28.kdams ,df_co_Syn_53.kdams ,df_co_Het_57.kdams , df_co_Het_58.kdams ,df_co_Het_59.kdams ,df_co_Het_60.kdams ] }, 
                                columns=['names','phis', 'kdams'])

fig4, (ax1)= plt.subplots(figsize = (7,7))
fig4.suptitle('Helper - coculture parameters', fontsize = 18)
ax1.set_xlabel('Detoxification rate ($\phi_{max})$',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$)',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

#ax1.set_xlim([1e-11, 1e-5]) #phi
#ax1.set_ylim([1e-5, 1e-2]) #kdam


for i, r in df_co_help.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = colors[i], marker = markers[i],markersize=11, linewidth=0.1, label=r['names'])

l4 = plt.legend(loc = 'best', fontsize = 12)

ax1.semilogy()
ax1.semilogx()

fig4.savefig('../figures/helper-tradeoffs')


############################################
#cPro params under coexistence
##############################################


colors = ['cornflowerblue','dodgerblue','steelblue','violet','blueviolet','mediumorchid','mediumpurple', 'lightgreen']
markers = ['o','o','o','v','v','v','v', 'D'] #HAVE THIS COORS[PMDWOTJ GENES]

df_co_pro = pd.DataFrame({'names' : ('Syn eWH7803_Vol52','Syn eWH7803_Vol28','Syn eCC9605_Vol53','Micromonas commoda - Vol_57','Micromonas pusilla - Vol_58','Ostreococcus lucimarinus - Vol_59','Ostreococcus tauri - Vol_60', 'Pro monoculture'),
                                'phis' : [df_co_Syn_52.phi,df_co_Syn_28.phi,df_co_Syn_53.phi,df_co_Het_57.phi, df_co_Het_58.phi,df_co_Het_59.phi,df_co_Het_60.phi, df_Pro_1.phi],
                                'kdams' : [ df_co_Syn_52.kdam ,df_co_Syn_28.kdam ,df_co_Syn_53.kdam ,df_co_Het_57.kdam , df_co_Het_58.kdam ,df_co_Het_59.kdam ,df_co_Het_60.kdam,df_Pro_1.kdam ] }, 
                                columns=['names','phis', 'kdams'])

fig5, (ax1)= plt.subplots(figsize = (7,7))
fig5.suptitle('Pro parameters in coculture', fontsize = 18)
ax1.set_xlabel('Detoxification rate ($\phi_{max})$',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$)',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

#ax1.set_xlim([1e-11, 1e-5]) #phi
#ax1.set_ylim([1e-5, 1e-2]) #kdam


for i, r in df_co_pro.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = colors[i], marker = markers[i],markersize=11, linewidth=0.1, label=r['names'])

l5 = plt.legend(loc = 'lower right', fontsize = 12)

ax1.semilogy()
ax1.semilogx()

fig5.savefig('../figures/Pro-helped')



print('Done!!!')

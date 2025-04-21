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


df_Pro_1 =  pd.read_csv('../data/inits/pro_MIT9215_inits4_1.csv') #Pro MIT9215
df_Pro_2 = pd.read_csv('../data/inits/pro_MIT9215_inits4_2.csv') #Pro MIT9215
df_Pros = pd.concat([df_Pro_1, df_Pro_2]).groupby(level=0).mean()   #avg of Pro MIT9215
df_Syn_28 = pd.read_csv('../data/inits/syn_vol28_inits4.csv') #Syn WH7803
df_Syn_52 = pd.read_csv('../data/inits/syn_vol52_inits4.csv') #Syn WH7803
df_Syn_WH7803 = pd.concat([df_Syn_28 , df_Syn_52]).groupby(level=0).mean()   #avg of Syn WH7803
df_Syn_53 = pd.read_csv('../data/inits/syn_vol53_inits4.csv') #Syn CC9605
df_Syn_54 = pd.read_csv('../data/inits/syn_vol54_inits4.csv') #Synechococcus WH8102 - no catalase strain  - growth doesn't occur in 0 H202
df_Het_57 = pd.read_csv('../data/inits/Het_57_inits4.csv') #Micromonas commoda
df_Het_58 = pd.read_csv('../data/inits/Het_58_inits4.csv') #Micromonas pusilla
df_Het_59 = pd.read_csv('../data/inits/Het_59_inits4.csv') #Ostreococcus lucimarinus
df_Het_60 = pd.read_csv('../data/inits/Het_60_inits4.csv') #Ostreococcus tauri
df_Het_55_3 = pd.read_csv('../data/inits/Het_55_inits_3-spike.csv') #Alteromonas macleodii  in 300nM 
df_Het_55_8 = pd.read_csv('../data/inits/Het_55_inits_8-spike.csv') #Alteromonas macleodii  in 800 nM
df_Het_55 = pd.concat([df_Het_55_3 , df_Het_55_8]).groupby(level=0).mean()   #avg of Alteromonas macleodii 

df_mono_no_ez55 = pd.DataFrame({'names' : ('Syn WH7803 avg' ,'Syn CC9605','Syn WH8102', 'Pro MIT9215 avg' ,'Micromonas commoda','Micromonas pusilla','Ostreococcus lucimarinus','Ostreococcus tauri'),
                                'phis' : [df_Syn_WH7803.phi,df_Syn_53.phi,df_Syn_54.phi,df_Pros.phi,df_Het_57.phi, df_Het_58.phi,df_Het_59.phi,df_Het_60.phi],
                                'kdams' : [df_Syn_WH7803.kdam,df_Syn_53.kdam,df_Syn_54.kdam,df_Pros.kdam,df_Het_57.kdam, df_Het_58.kdam,df_Het_59.kdam,df_Het_60.kdam ], 
                                'markers' : ['o','o','o','D','v','v','v','v'],
                                'colors' : ['dodgerblue','steelblue','b','lightgreen','violet','blueviolet','mediumorchid','mediumpurple']}, 
                                columns=['names','phis', 'kdams','markers','colors'])

df_mono = pd.DataFrame({'names' : ('$Synechococcus$ WH7803' ,'$Synechococcus$ CC9605','$Synechococcus$ WH8102', '$Prochlorococcus$ MIT9215' ,'$Micromonas$ $commoda$','$Micromonas$ $pusilla$','$Ostreococcus$ $lucimarinus$','$Ostreococcus$ $tauri$', '$Alteromonas$ $macleodii$'),
                                'phis' : [df_Syn_WH7803.phi,df_Syn_53.phi,df_Syn_54.phi,df_Pros.phi,df_Het_57.phi, df_Het_58.phi,df_Het_59.phi,df_Het_60.phi,df_Het_55.phi],
                                'kdams' : [df_Syn_WH7803.kdam,df_Syn_53.kdam,df_Syn_54.kdam,df_Pros.kdam,df_Het_57.kdam, df_Het_58.kdam,df_Het_59.kdam,df_Het_60.kdam ,df_Het_55.kdam], 
                                'markers' : ['o','o','o','D','v','v','v','v','d'],
                                'colors' : ['dodgerblue','steelblue','b','lightgreen','violet','blueviolet','mediumorchid','mediumpurple','k']}, 
                                columns=['names','phis', 'kdams','markers','colors'])



df_co_Syn_28 = pd.read_csv('../data/inits/coculture_1,28_pro_inits4.csv')
df_co_Syn_52 = pd.read_csv('../data/inits/coculture_1,52_pro_inits4.csv')
df_co_Syn_WH7803 = pd.concat([df_co_Syn_28 , df_co_Syn_52]).groupby(level=0).mean()   #avg of Syn WH7803

df_co_Syn_53 = pd.read_csv('../data/inits/coculture_1,53_pro_inits4.csv') #Syn CC9605 in coculture with Pro MIT9215
df_co_Syn_54 = pd.read_csv('../data/inits/coculture_1,54_pro_inits4.csv') #Syn WH8102 in coculture with Pro MIT9215

df_co_Het_57 = pd.read_csv('../data/inits/Het_coculture_1,57_pro_inits4.csv')
df_co_Het_58 = pd.read_csv('../data/inits/Het_coculture_1,58_pro_inits4.csv')
df_co_Het_59 = pd.read_csv('../data/inits/Het_coculture_1,59_pro_inits4.csv')
df_co_Het_60 = pd.read_csv('../data/inits/Het_coculture_1,60_pro_inits4.csv')


df_co_help = pd.DataFrame({'names' : ('Syn WH7803_avg','Syn CC9605 ','Syn WH8102' ,'Micromonas commoda','Micromonas pusilla','Ostreococcus lucimarinus','Ostreococcus tauri'),
                                'phis' : [df_co_Syn_WH7803.phis,df_co_Syn_53.phis,df_co_Syn_54.phis,df_co_Het_57.phis, df_co_Het_58.phis,df_co_Het_59.phis,df_co_Het_60.phis],
                                'kdams' : [ df_co_Syn_WH7803.kdams ,df_co_Syn_53.kdams ,df_co_Syn_54.kdams,df_co_Het_57.kdams , df_co_Het_58.kdams ,df_co_Het_59.kdams ,df_co_Het_60.kdams ], 
                                'markers' : ['o','o','o','v','v','v','v'],
                                'colors' : ['dodgerblue','steelblue','b','violet','blueviolet','mediumorchid','mediumpurple']}, 
                                columns=['names','phis', 'kdams','markers','colors'])


df_co_pro = pd.DataFrame({'names' : ('Syn WH7803 avg','Syn CC9605 ','Syn WH8102','Micromonas commoda','Micromonas pusilla','Ostreococcus lucimarinus','Ostreococcus tauri', 'Pro monoculture'),
                                'phis' : [df_co_Syn_WH7803.phi,df_co_Syn_53.phi,df_co_Syn_54.phi,df_co_Het_57.phi, df_co_Het_58.phi,df_co_Het_59.phi,df_co_Het_60.phi, df_Pros.phi],
                                'kdams' : [ df_co_Syn_WH7803.kdam ,df_co_Syn_53.kdam ,df_co_Syn_54.kdam,df_co_Het_57.kdam , df_co_Het_58.kdam ,df_co_Het_59.kdam ,df_co_Het_60.kdam,df_Pros.kdam ] , 
                                'markers' : ['o','o','o','v','v','v','v','D'],
                                'colors' : ['dodgerblue','steelblue','b','violet','blueviolet','mediumorchid','mediumpurple','lightgreen']}, 
                                columns=['names','phis', 'kdams','markers','colors'])

####################################################
#monoculture Pro, Syn only
####################################################


df_mini = pd.DataFrame({'names' : ('$Synechococcus$ WH7803 ' ,'$Synechococcus$ CC9605', '$Prochlorococcus$ MIT9215' ,'$Alteromonas$ $macleodii$'),
                                'phis' : [df_Syn_WH7803.phi,df_Syn_53.phi,df_Pros.phi,df_Het_55.phi],
                                'kdams' : [df_Syn_WH7803.kdam,df_Syn_53.kdam,df_Pros.kdam,df_Het_55.kdam],
                                'markers' : ['o','o','D','d'],
                                'colors' : ['dodgerblue','steelblue','lightgreen','k']}, 
                                columns=['names','phis', 'kdams','markers','colors'])


##mini fig of monoculture tradeoffs  - Manuscript of CHap 2
fig3, (ax1)= plt.subplots(figsize = (7,6))
ax1.set_xlabel('Detoxification rate ($\phi_{det,i})$',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$)',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

#ax1.set_xlim([1e-9, 1e-5]) #phi
#ax1.set_ylim([1e-5, 1e-2]) #kdam


for i, r in df_mini.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = r['colors'], marker = r['markers'],markersize=11, linewidth=0.1, label=r['names'])

l3 = plt.legend(loc = 'best', fontsize = 12)

ax1.semilogy()
ax1.semilogx()

#fig3.savefig('../figures/tradeoffs_mini')


####################################################
#monoculture of all  Pro, Syn, and Het Kdam and Phis
#####################################################

#Not all of the data sets have cocurrent HOOH measurements so Phi is less constrained and therefore can be biassing the conclusions a lot. 
fig3, (ax1)= plt.subplots(figsize = (7,6))
ax1.set_xlabel('Detoxification rate ($\phi_{det,i})$',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam,i}$)',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)


for i, r in df_mono.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = r['colors'], marker = r['markers'],markersize=11, linewidth=0.1, label=r['names'])

l3 = plt.legend(loc = 'best', fontsize = 12)

l3.draw_frame(False)
ax1.semilogy()
ax1.semilogx()

fig3.subplots_adjust(bottom=0.15)

fig3.savefig('../figures/tradeoffs')



############################################
#coexistence Phis and Kdams
##############################################


fig4, (ax1)= plt.subplots(figsize = (7,6))
fig4.suptitle('Helper - coculture parameters', fontsize = 18)
ax1.set_xlabel('Detoxification rate ($\phi_{det,i})$ of coexistant microbe',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$) of coexistant microbe',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

#ax1.set_xlim([1e-11, 1e-5]) #phi
#ax1.set_ylim([1e-5, 1e-2]) #kdam


for i, r in df_co_help.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = r['colors'], marker = r['markers'],markersize=11, linewidth=0.1, label=r['names'])

l4 = plt.legend(loc = 'best', fontsize = 12)



ax1.semilogy()
ax1.semilogx()

#fig4.savefig('../figures/helper-tradeoffs')


############################################
#cPro params under coexistence
##############################################



fig5, (ax1)= plt.subplots(figsize = (7,6))
fig5.suptitle('Pro parameters in coculture', fontsize = 18)
ax1.set_xlabel('Detoxification rate ($\phi_{det,i})$ of Pro ',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$) of Pro ',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

#ax1.set_xlim([1e-11, 1e-5]) #phi
#ax1.set_ylim([1e-5, 1e-2]) #kdam


for i, r in df_co_pro.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = r['colors'], marker = r['markers'],markersize=11, linewidth=0.1, label=r['names'])

l5 = plt.legend(loc = 'lower right', fontsize = 12)

ax1.semilogy()
ax1.semilogx()

#fig5.savefig('../figures/Pro-helped')





####################
#ratios
#####################

r_28 = (df_Syn_28.phi)/(df_co_Syn_28.phis)
r_52 = (df_Syn_52.phi)/(df_co_Syn_52.phis)
r_53 = (df_Syn_53.phi)/(df_co_Syn_53.phis)


r_57 = (df_Het_57.phi)/(df_co_Het_57.phis)
r_58 = (df_Het_58.phi)/(df_co_Het_58.phis)
r_59 = (df_Het_59.phi)/(df_co_Het_59.phis)
r_60 = (df_Het_60.phi)/(df_co_Het_60.phis)


phi_rs = [r_28,r_52,r_53,r_57,r_58,r_59,r_60]

P = df_Pro_1.kdam

r_28 = (P)/(df_co_Syn_28.kdam)
r_52 = (P)/(df_co_Syn_52.kdam)
r_53 = (P)/(df_co_Syn_53.kdam)


r_57 = (P)/(df_co_Het_57.kdam)
r_58 = (P)/(df_co_Het_58.kdam)
r_59 = (P)/(df_co_Het_59.kdam)
r_60 = (P)/(df_co_Het_60.kdam)


kdamp_rs = [r_28,r_52,r_53,r_57,r_58,r_59,r_60]

colors = ['dodgerblue','cornflowerblue','steelblue','violet','blueviolet','mediumorchid','mediumpurple']
markers = ['o','o','o','v','v','v','v'] #HAVE THIS COORS[PMDWOTJ GENES]


df_ratios = pd.DataFrame({'names' : ('Syn WH7803_Vol28','Syn WH7803_Vol52','Syn CC9605 ','Micromonas commoda','Micromonas pusilla','Ostreococcus lucimarinus','Ostreococcus tauri'),
                                'phis' :phi_rs ,
                                'kdams' : kdamp_rs }, 
                                columns=['names','phis', 'kdams'])



fig6, (ax1)= plt.subplots(figsize = (7,6))
fig6.suptitle('Helper ratios', fontsize = 18)
ax1.set_xlabel('((P_mono kdam) /(P_co kdam))',fontsize = 16)
ax1.set_ylabel('((Helper_mono phi)/(Helper_co phi))',fontsize = 16)

for i, r in df_ratios.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = colors[i], marker = markers[i],markersize=11, linewidth=0.1, label=r['names'])

l6 = plt.legend(loc = 'lower right', fontsize = 12)

ax1.semilogy()
ax1.semilogx()

#fig6.savefig('../figures/hlper_ratios')



print('Done!!!')

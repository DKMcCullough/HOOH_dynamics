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

df = pd.DataFrame({'names' : ('Syn eWH7803_Vol52','Syn eWH7803_Vol28','Syn eCC9605_Vol53','Pro eMIT9215_Vol1_sample1','Pro eMIT9215_Vol1_sample2'),
                                'phis' : [df_Syn_52.phi,df_Syn_28.phi,df_Syn_53.phi,df_Pro_1.phi,df_Pro_2.phi],
                                'kdams' : [df_Syn_52.kdam,df_Syn_28.kdam,df_Syn_53.kdam,df_Pro_1.kdam,df_Pro_2.kdam] }, 
                                columns=['names','phis', 'kdams'])


#Not all of the data sets have cocurrent HOOH measurements so Phi is less constrained and therefore can be biassing the conclusions a lot. 
fig3, (ax1)= plt.subplots(figsize = (7,6))
fig3.suptitle('Tradeoffs', fontsize = 17)
ax1.set_xlabel('Detoxificatiion rate ($\phi_{max})$',fontsize = 14)
ax1.set_ylabel('Damage rate ($\kappa_{dam}$)',fontsize = 14)
ax1.tick_params(axis = 'both', which = 'both', length = 4, labelsize = 12)


ax1.set_xlim([0.000000000001, 0.0001]) #phi
ax1.set_ylim([0.0003, 0.009]) #kdam


colors = ['cornflowerblue','dodgerblue','steelblue','g','lightgreen']
markers = ['o','o','o','D','D'] #HAVE THIS COORS[PMDWOTJ GENES]

#df.plot.scatter(x='kdams', y='phis',xlabel = 'kdam', ylabel = 'phi',label = 'names')
for i, r in df.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = colors[i], marker = markers[i],markersize=10, linewidth=0.1, label=r['names'])

#plt.axes(fontsize = 12)
l3 = plt.legend(loc = 'lower left', fontsize = 10)
#l3.draw_frame(False)




ax1.semilogy()
ax1.semilogx()

fig3.savefig('../figures/tradeoffs')


print('Done!!!')

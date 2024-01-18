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



df = pd.DataFrame({'names' : ('Syn WH7803','Pro MIT9215'),
                                'phis' : [1.71E-06,2.82E-12],
                                'kdams' : [0.00132393992680983,0.00322116921322925] }, 
                                columns=['names','phis', 'kdams'])

#Not all of the data sets have cocurrent HOOH measurements so Phi is less constrained and therefore can be biassing the conclusions a lot. 
fig3, (ax1)= plt.subplots(figsize = (7,5))
fig3.suptitle('Tradeoffs', fontsize = 17)
ax1.set_ylabel('phi',fontsize = 15)
ax1.set_xlabel('kdam',fontsize = 15)
ax1.tick_params(axis = 'both', which = 'both', length = 6, labelsize = 14)


#ax1.set_xticklabels(labels = 'kdam' ,fontsize = 14)
#ax1.set_yticklabels(labels = 'phi' ,fontsize = 14)

colors = ['orange','g']
#shape = [+,-,+,-,-,-,-] #HAVE THIS COORS[PMDWOTJ GENES]

#df.plot.scatter(x='kdams', y='phis',xlabel = 'kdam', ylabel = 'phi',label = 'names')
for i, r in df.iterrows():
    ax1.plot(r['kdams'], r['phis'], 'o',c = colors[i], markersize=7, linewidth=0.1, label=r['names'])

#plt.axes(fontsize = 12)
plt.legend(loc='best')



ax1.semilogy()
ax1.semilogx()

plt.show()

fig3.savefig('../figures/tradeoffs')


print('Done!!!')
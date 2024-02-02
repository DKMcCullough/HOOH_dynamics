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

df = pd.DataFrame({'names' : ('Syn WH7803_vol52','Syn WH7803_vol28','Syn CC9605_vol53','Pro MIT9215_sample1','Pro MIT9215_sample2'),
                                'phis' : [1.94E-06,1.64E-06,1.61E-06,3.95E-15,1.87E-15],
                                'kdams' : [0.001390782,0.001691207,0.000656184,0.003121548,0.005126291] }, 
                                columns=['names','phis', 'kdams'])

#Not all of the data sets have cocurrent HOOH measurements so Phi is less constrained and therefore can be biassing the conclusions a lot. 
fig3, (ax1)= plt.subplots(figsize = (8,6))
fig3.suptitle('Tradeoffs', fontsize = 17)
ax1.set_xlabel('\u03C6',fontsize = 14)
ax1.set_ylabel('kdam',fontsize = 14)
ax1.tick_params(axis = 'both', which = 'both', length = 4, labelsize = 12)


ax1.set_xlim([0.0000000000000001, 0.0001]) #phi
ax1.set_ylim([0.0001, 0.01]) #kdam


colors = ['cornflowerblue','dodgerblue','steelblue','g','lightgreen']
markers = ['o','o','o','d','d'] #HAVE THIS COORS[PMDWOTJ GENES]

#df.plot.scatter(x='kdams', y='phis',xlabel = 'kdam', ylabel = 'phi',label = 'names')
for i, r in df.iterrows():
    ax1.plot(r['phis'], r['kdams'], 'o',c = colors[i], marker = markers[i],markersize=10, linewidth=0.1, label=r['names'])

#plt.axes(fontsize = 12)
l3 = plt.legend(loc = 'lower left', fontsize = 10)
#l3.draw_frame(False)




ax1.semilogy()
ax1.semilogx()

plt.show()

fig3.savefig('../figures/tradeoffs')


print('Done!!!')
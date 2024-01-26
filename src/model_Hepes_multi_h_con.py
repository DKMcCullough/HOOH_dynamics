'''
All Hepes trials from Fig 1_a

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Goal is to have all Hepes data shown and then have analytical solutions for all that differ only in production (S_hooh)terms...with decay staying the same for all.

Using Morrisand Zinser 2013 Fig 1_a  data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_5_Hepes_Morris_2013/scripts/
'''


from scipy import *
from scipy.integrate import odeint
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   
import scipy 
import ODElib
import random as rd
import sys



######################################################################

#  Importing the data as dfs for each trial

######################################################################
df_all = pd.read_csv("../data/Hepes_assay_2013.csv") #use relative paths 
df_all.rename(columns = {'time(days)':'time'},inplace = True) 
df_all.dropna()

df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])
df_all['log4'] = np.log(df_all['rep4'])
df_all['log5'] = np.log(df_all['rep5'])
df_all['log6'] = np.log(df_all['rep6'])

df_all['avg1'] = np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['avg2'] = np.nanmean(np.r_[[df_all[i] for i in ['rep4', 'rep5', 'rep6']]],axis=0)
df_all['stdv1'] = np.std(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['stdv2'] = np.std(np.r_[[df_all[i] for i in ['rep4', 'rep5', 'rep6']]],axis=0)

#bio rep logged avgs and stdvs 
df_all['lavg1'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['lavg2'] = np.nanmean(np.r_[[df_all[i] for i in ['log4','log5','log6']]],axis=0)
df_all['stdlog1'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['stdlog2'] = np.std(np.r_[[df_all[i] for i in ['log4','log5','log6']]],axis=0)




#calculating stats for trial compilation
df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4','rep5','rep6']]],axis=0)
df_all['sigma'] = np.nanstd(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4','rep5', 'rep6']]],axis=0)

df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3', 'log4','log5','log6']]],axis=0)
df_all['log_sigma'] = np.nanstd(np.r_[[df_all[i] for i in ['log1','log2','log3','log4','log5','log6']]],axis=0)




ROSs = df_all['Hepes_treatment'].unique()
nROSs = ROSs.shape[0]
#ROSs.sort()     #Change here to get ambient ROS to be 61 or 90 for species and 0 is detoxed. 


f1,ax1 = plt.subplots(figsize=[8,6])

colors = ('yellowgreen','green','c','b','purple','orange','r','k')
markers = ('.','P','s','1','^','*','p','+')


h0s = np.array([])

for (ros,ri) in zip(ROSs,range(nROSs)): # loop over ROS
    tdf = (df_all[df_all['Hepes_treatment']==ros]) # select one ROS treatment at a time 
    hoohs = tdf['log_abundance']
    h0 = hoohs.iloc[0]
    h0s = np.append(h0,h0s) #fing start H for each trial for easier modeling later
    ax1.plot(tdf['time'],tdf['log_abundance'], marker='s',color = colors[ri], markersize = 10, linestyle = ':', label='HEPES = '+str(ros)) # graph each 







plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel('HOOH concentration (\u03BCM)', fontsize = 20)
plt.yscale('log')
plt.tick_params(labelsize = 18) 
plt.tick_params(size = 14)
plt.legend(loc='upper left')
#plt.show()



####################################

 # Set Up model parts 

####################################
    
step = 0.05 #delta t
ndays = 26
times = np.linspace(0,ndays,int(ndays/step))

h_cons = np.array([0.15,0.6,1.1,1.1])
#h_convert = 0.6     #term used to convert input Hepes concentration to HOOH felt by cells
S_Hs = np.array([])
'''
#making HOOH [ ] fromo HEPES [ ] given 
for r in ROSs: 
    S_HOOH = r * h_convert      #TAKING HEPES AT EACH TREATMENT AND FINDING HOOH BASED ON COMMON CONVERSION FACTOR
    S_Hs = np.append(S_Hs,S_HOOH)
'''
for (h_con,ros) in zip(h_cons,ROSs): # loop over ROS
    S_HOOH = ros * h_con
    print(ros,S_HOOH)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    S_Hs = np.append(S_Hs,S_HOOH)

nS_Hs = S_Hs.shape[0]


#model for use in ode int
def HsODEint(y,t,S_HOOH,delta):
    H = y[0]
    dHdt = S_HOOH-delta*H
    #print(dHdt,t,H)
    return dHdt 

liROS = ROSs.tolist()

delta = 0.02
'''
for (ros,S_HOOH,h0) in zip(ROSs,S_Hs,h0s): # loop over ROS
    count = (liROS.index(ros))
    print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,h0,times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')
'''
for (ros,S_HOOH) in zip(ROSs,S_Hs): # loop over ROS
    count = (liROS.index(ros))
    #print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,r_[[0.7]],times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')

plt.xlim([0, 26])
plt.ylim([0.15, 50])

'''
####################################

#analytical solutions 

####################################


S_HOOH = S_Hs

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 

plt.plot(times,Hs,c='y',marker='.',label='Analytical Solution')



delta = 0.17
for (ros,S_HOOH) in zip(ROSs,S_Hs): # loop over ROS
    count = (liROS.index(ros))
    #print(count)
    solutions = Hs = f(times,S_HOOH,delta)
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')

'''

plt.legend()
plt.show()

fig, (ax2,ax3) = plt.subplots(1,2)
fig.suptitle('HEPES conversions',fontsize = 20)
plt.subplots_adjust(wspace = 0.3, top = 0.85)

#f2, (ax2, ax3) = plt.subplots(1,2, figsize=[8,6])
ax2.plot(ROSs,S_Hs, marker = 's', markersize = 8, color = 'brown')
ax2.set_title('HEPES vs modeled S_HOOH',fontsize = 10)
ax2.set_xlabel('HEPES added', fontsize = 10)
ax2.set_ylabel('Modeled HOOH Supply', fontsize = 10)

ax3.plot(ROSs,h_cons, marker = 'd', markersize = 8, color = 'k')
ax3.set_title('h_convert',fontsize = 10)
ax3.set_xlabel('HEPES added', fontsize = 10)
ax3.set_ylabel('h_convert needed for model', fontsize = 10)
#plt.tick_params(labelsize = 18) 
#plt.tick_params(size = 14)
#plt.legend(loc='lower right')






# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')
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




######################################################################

#  Importing the data as dfs for each trial

######################################################################

import csv
df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'hepes_array-2013Morrisfig2a', header  = 0)

#df_all = pd.read_csv("../data/Hepes_assay_2013.csv") #use relative paths 
df_all.rename(columns = {'time(days)':'time'},inplace = True) 
df_all.dropna()
#use df clip? or just use 0H as control for base level of media 

#making log of data to look at error (how to do in loop? all combined whe I tried)
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


#split ROS treatments by number 
ROSs = df_all['Hepes_treatment'].unique()
ROSs = ROSs[~pd.isna(ROSs)] #clipping off nans 
nROSs = ROSs.shape[0]


#model set up 
    
step = 0.05 #delta t
ndays = 26
times = np.linspace(0,ndays,int(ndays/step))

#h_cons = np.array([0.15,0.6,1.1,1.1])
h_convert = 60     #term used to convert input Hepes concentration to HOOH felt by cells
S_Hs = np.array([])

#making HOOH [ ] fromo HEPES [ ] given 
for r in ROSs: 
    S_HOOH = r * h_convert      #TAKING HEPES AT EACH TREATMENT AND FINDING HOOH BASED ON COMMON CONVERSION FACTOR
    S_Hs = np.append(S_Hs,S_HOOH)


nS_Hs = S_Hs.shape[0]
deltas = np.array([0.4,0.4,0.04,0.04,0.04,0.0,0.0])

#model for use in ode int
def HsODEint(y,t,S_HOOH,delta):
    H = y[0]
    dHdt = S_HOOH-delta*H
    #print(dHdt,t,H)
    return dHdt 

liROS = ROSs.tolist()
h0s = np.array([])

#lists for loop 
colors = ('yellowgreen','green','c','b','purple','orange','r','k')
ls = ('solid','dotted','dashed','dashdot','solid','dotted','dashed','dashdot')


f1,ax1 = plt.subplots(figsize=[10,5])

for (ros,ri,S_HOOH,delta) in zip(ROSs,range(nROSs),S_Hs,deltas): # loop over ROS
    tdf = (df_all[df_all['Hepes_treatment']==ros]) # select one ROS treatment at a time 
    hoohs = tdf.abundance
    h0 = hoohs.iloc[0]
    h0s = np.append(h0,h0s) #fing start H for each trial for easier modeling later
    ax1.plot(tdf['time'],tdf['abundance'], marker = 's',color = colors[ri], markersize = 10, linestyle = ':', label='HEPES (\u03BCM) = '+str(ros))
    ax1.errorbar(tdf['time'],tdf['abundance'],yerr = tdf['sigma'],marker = 's',color = colors[ri])
    count = (liROS.index(ros))
    #print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,(np.full(1,7)),times,args=(S_HOOH,delta))  
    ax1.plot(times,solutions[:,0],color = colors[count],linestyle = ls[count], marker = 'None')



'''
for (ros,S_HOOH,delta) in zip(ROSs,S_Hs,deltas): # loop over ROS
    count = (liROS.index(ros))
    #print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,(np.full(1,7)),times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count],linestyle = ls[count], marker = 'None')
'''





plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel('HOOH concentration (\u03BCM)', fontsize = 20)
plt.yscale('log')
plt.tick_params(labelsize = 18) 
plt.tick_params(size = 14)
#plt.legend(loc='upper left')
#plt.show()



####################################

 # Set Up model parts 

####################################




plt.xlim([0, 26])
plt.ylim([0.15, 50])



plt.legend()
plt.show()


fig, (ax3) = plt.subplots()
fig.suptitle('Delta dynamics',fontsize = 20)
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax3.plot(ROSs,deltas, marker = 'd', markersize = 8, color = 'k')
ax3.set_title('delta',fontsize = 10)
ax3.set_xlabel('HEPES added', fontsize = 10)
ax3.set_ylabel('delta needed for model', fontsize = 10)
#plt.tick_params(labelsize = 18) 
#plt.tick_params(size = 14)
#plt.legend(loc='lower right')






# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')
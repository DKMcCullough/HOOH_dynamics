'''

name:   ma_2018_all_graphed.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Project_3_Temperature_interactions/src/ma_2018_all_graph.py'

author: DKM


goal: import and all treatments of NH4 addition experients using MIT9215


working on: Get params into loops so you can have differfent kdams for differentt lines


'''



import pandas as pd
import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
from scipy.integrate import *
from scipy import *
from pylab import *



############################

#  Data Import from csv   

############################

df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'Ma_2018', header = 0)

#df_all = pd.read_csv("../data/Ma_fig_1_raw_data.csv") #use relative paths 
df_all[['rep1', 'rep2','rep3']] = df_all[['rep1', 'rep2','rep3']].fillna(value=0)
df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first
df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'
#neeed to cut off extra NAN columns 

master_df = df_all 

strains = master_df['strain'].unique()
temps = master_df['temperature'].unique()
my_temps = np.r_[[temps[i] for i in [0,1,3,4,5,6,7]]]
ROSs = master_df['treatment'].unique()
ROSs.sort()     #Change here to get ambient ROS to be 61 or 90 for species and 0 is detoxed. 


####################################

# Slicing Data

####################################

'''
#df_all[df_all['strain'] == 64]

dc = dict()
####slicing df into treatments and saving into dictionary (dc)######
for i in bugs:
    df_i = df_all[df_all["strain"].isin([i])]    
    for i in treatments:
        df_i = df_all[df_all["treatment"].isin([i])]
        name = ('df_' + str(i))
        dc.update({name : df_i})  #update dictionary of dfs with each loop itteration. 


avgs = dict()
yerrs = dict()

#calculating avg point and yerr for each treatment

for i in dc :
    df_i = dc[i]   
    rep_df_i = df_i[['rep1', 'rep2', 'rep3']]
    avg_i = rep_df_i.mean(axis=1)  
    avgs.update({'avg_'+i : avg_i })
    yerr_i = rep_df_i.std(axis=1)   
    yerrs.update({'yerr_'+i : yerr_i })


####################################

#graphing 

####################################
'# read in data
master_df = pd.read_csv('../data/ma_temp_ros.csv')
master_df['abundance'] = np.nanmean(np.r_[[master_df[i] for i in ['rep1','rep2','rep3']]],axis=0)
'''
# unique strings
strains = master_df['strain'].unique()
temps = master_df['temperature'].unique()
my_temps = np.r_[[temps[i] for i in [0,1,3,4,5,6,7]]]
ROSs = master_df['treatment'].unique()
ROSs.sort()     #Change here to get ambient ROS to be 61 or 90 for species and 0 is detoxed. 
master_df['abundance'] = np.nanmean(np.r_[[master_df[i] for i in ['rep1','rep2','rep3']]],axis=0)
#my_ROSs = (0.0,ambient (61 or 90)),200,400 

# number of treatments
nstrains = strains.shape[0]
ntemps = my_temps.shape[0]
nROSs = ROSs.shape[0]

# setup figures
f1,ax1 = plt.subplots(ntemps,nstrains,figsize=[6,12])
f2,ax2 = plt.subplots(nROSs,nstrains,figsize=[6,12])
f3,ax3 = plt.subplots(ntemps,nROSs,figsize=[12,12])

# plot
for (strain,si) in zip(strains,range(nstrains)): # first loop over strains
    for (temp,ti) in zip(my_temps,range(ntemps)): # now loop over temps
        for (ros,ri) in zip(ROSs,range(nROSs)): # now loop over ROS
            tdf = master_df[(master_df['strain']==strain) \
                & (master_df['temperature']==temp) \
                & (master_df['treatment']==ros)] # select dataset
            ax1[ti,si].plot(tdf['times'],tdf['abundance'],marker='o',label='ros ='+str(ros)) # ROS varies, temp constant
            ax2[ri,si].plot(tdf['times'],tdf['abundance'],marker='o',label='temp ='+str(temp)) # ROS constant, temp varies
            ax3[ti,ri].plot(tdf['times'],tdf['abundance'],marker='o',label=strain) # ROS and temp vary
            ax1[ti,si].semilogy() # logscale
            ax2[ri,si].semilogy()
            ax3[ti,ri].semilogy()
            if (si == nstrains-1): # annotate
                ax1[ti,1].text(1.2,0.6,'Temp: '+ '\n' + str(temp)+' ($^\circ$C)',horizontalalignment='center',\
                        verticalalignment='center', transform=ax1[ti,1].transAxes)
                ax2[ri,1].text(1.2,0.6,'ROS: '+'\n' + str(ros)+' nM',horizontalalignment='center',\
                        verticalalignment='center', transform=ax2[ri,1].transAxes)
                ax3[ti,-1].text(1.2,0.4,'Temp: '+ '\n' +str(temp)+' ($^\circ$C)',horizontalalignment='center',\
                        verticalalignment='center', transform=ax3[ti,-1].transAxes)
                ax3[0,ri].set_title('ROS: '+str(ros)+' nM')

# legends
l1 = ax1[0,0].legend(ncol=2)
l2 = ax2[1,0].legend()
l3 = ax3[0,0].legend()
l1.draw_frame(False)
l2.draw_frame(False)
l3.draw_frame(False)

# make space on the right for annotation (e.g. ROS=0, etc.)
f1.subplots_adjust(right=0.85, wspace = 0.3, hspace = 0.45)
f2.subplots_adjust(right=0.85, wspace = 0.3, hspace = 0.45)
f3.subplots_adjust(right=0.9, wspace = 0.3, hspace = 0.45)

# titles
ax1[0,0].set_title(strains[0])
ax1[0,1].set_title(strains[1])
ax2[0,0].set_title(strains[0])
ax2[0,1].set_title(strains[1])

# xlabels
for a in ax1[-1,:]:
    a.set_xlabel('Time (days)')
for a in ax2[-1,:]:
    a.set_xlabel('Time (days)')
for a in ax3[-1,:]:
    a.set_xlabel('Time (days)')

# ylabels
for a in ax1[:,0]:
    a.set_ylabel('Cells (ml$^{-1}$)')
for a in ax2[:,0]:
    a.set_ylabel('Cells (ml$^{-1}$)')
for a in ax3[:,0]:
    a.set_ylabel('Cells (ml$^{-1}$)')

f1.savefig('../figures/temps')
f2.savefig('../figures/ross')
f3.savefig('../figures/temp_x_ros')


### To do 

#detoxed, ambient, 200, 400 

'''
my_ROSs = np.r_[[]]
for r in ROSs: 
    if r < 50:
        ROS = 'detoxed'
    if 50 < r < 100:
        ROS = 'ambient'
    if r>100: 
        ROS = r
    my_ROSs = np.append(my_ROSs, ROS)


u_ROSs,indices = np.unique(my_ROSs, return_index = True) #gets indicies of unique values (getting rid of extra ambient)
#print(u_ROSs)
#print(indices)
indices.sort()
#print(indices)
my_ROSs = (my_ROSs[indices])
print(my_ROSs)
'''
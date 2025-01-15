'''

name:   viz_spikes.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser/Projects/ROS.../HOOH.....src'

author: DKM



'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import seaborn as sns

#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'

######################################################
#reading in data and configureing 
#####################################################
df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_1-31-dataset', header = 1)
df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

# hydrogen peroxide abiotic spikes
df_a = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()
df0_abiotic = df_a.loc[~ df_a['assay'].str.contains('4', case=False)]  #assay 0 H
df4_abiotic = df_a.loc[(df_a['assay'].str.contains('4', case=False))]

# syn
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()

# pro with h202
df0_pro = df_mono.loc[~ df_mono['assay'].str.contains('4', case=False) & (df_mono['Vol_number']== 1)]  #assay 0 H 
df4_pro = df_mono.loc[(df_mono['assay'].str.contains('4', case=False)) & (df_mono['Vol_number']== 1)]

# syn with h202
df0_syn = df_mono.loc[~ df_mono['assay'].str.contains('4', case=False) & (df_mono['Vol_number']== 52)]  #assay 0 H
df4_syn = df_mono.loc[(df_mono['assay'].str.contains('4', case=False)) & (df_mono['Vol_number']== 52)]






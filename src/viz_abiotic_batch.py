'''

name:   viz_abiotic_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model of Monoculture BCC assays to graph 0 H phyotplankton biomass 

working on: - getting this in model to play so we can model all at once 

'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 


#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'

######################################################
#reading in data and configureing 
#####################################################

df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
#splicing abiotic and mono or coculture data
df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  

#setting working directory 
#only reps 1 and 3
df = df_abiotic

df['log1'] = np.log(df['rep1'])
df['log3'] = np.log(df['rep3'])


#total avgs and stdvs
df['abundance'] =  np.nanmean(np.r_[[df[i] for i in ['rep1','rep3']]],axis=0)
df['sigma'] = np.std(np.r_[[df[i] for i in ['rep1','rep3']]],axis=0)
df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['log1','log3']]],axis=0)
df['log_sigma'] = np.std(np.r_[[df[i] for i in ['log1','log3']]],axis=0)



#picking 0 or 400 assay dfs
#df0 = df.loc[df['assay'].str.contains('_0', case=False)].copy()
#df400 = df.loc[df['assay'].str.contains('_0', case=False)].copy()

treats = df['assay'].unique()
ntreats = treats.shape[0]


    #setting up fig 1 for dyanmics 


#####################################################
# Set up large loop  (treatments) 
#####################################################


fig1,ax1 = plt.subplots(2,2,figsize=[13,9]) #plot creation and config 
fig1.suptitle('Raw dynamics') #full title config
ax1[0,0].set_title('HOOH Concentration')
ax1[0,1].set_title('avg vs stdv')
fig1.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.85, wspace=0.25, hspace=0.25) #shift white space for better fig view
ax1[0,0].set_xlabel('Time (days)')
ax1[0,0].set_ylabel('HOOH concentration (\u03BCM)')
ax1[0,1].set_xlabel('Raw Mean')
ax1[0,1].set_ylabel('Raw standard deviation ')

fig2,ax2 = plt.subplots(2,2,figsize=[12,10]) #plot creation and config 
#full title config
fig2.suptitle('Logged dynamics ')
ax2[0,0].set_title('HOOH concentration ')
ax2[0,1].set_title('Log avg vs stdv')
fig2.subplots_adjust(left=0.1, bottom=0.10, right=0.85, top=0.85, wspace=0.25, hspace=0.35) #shift white space for better fig view
ax2[0,0].set_xlabel('Time (days)')
ax2[0,0].set_ylabel('HOOH concentration (\u03BCM)')
ax2[0,1].set_xlabel('Log Mean')
ax2[0,1].set_ylabel('Log Standard deviation')

fig3,ax3 = plt.subplots(2,2,figsize = [8,6])
fig3.suptitle("Pearsons's R coorelations")
fig3.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.85, wspace=0.25, hspace=0.4)
ax3[0,0].set_title('( 0 abiotic) ')
ax3[0,1].set_title('( 400 abiotic) ')
    
#####################################################
# Set up loop of vol numbers inside treatment loop  
#####################################################

for (t,nt) in zip(treats,range(ntreats)):

    df = df_abiotic[((df_abiotic['assay']==t))].copy()
    #ratios for graphing stretch
    raw_ratio =  (np.max(df.abundance/df.sigma))
    log_ratio =  (np.max(df.log_abundance/df.log_sigma))
    rats  = (raw_ratio,log_ratio)
    rats = np.rint(rats)
    best_ratio = (np.max(rats))
    print(best_ratio)
    rat_dif = (rats/best_ratio)
    stretch_rat = abs(rat_dif -1)
    print(stretch_rat)
#setting working df as a single Experiment in df_all
    ax1[nt,0].plot(df.time,df.rep1,color = 'pink', label = 'rep1')
    ax1[nt,0].plot(df.time,df.rep3, color = 'yellow',label = 'rep3')
    ax1[nt,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker = '*', c='purple',label =  'Mean')
    ax1[nt,1].scatter(df.abundance,df.sigma, c='purple')
    ax1[nt,1].text(1.2,0.5,str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
    ax1[nt,0].semilogy()

    ax2[nt,0].plot(df.time,df.log1, color = 'b', label = 'l1')
    ax2[nt,0].plot(df.time,df.log3, color = 'g',label = 'l3')
    ax2[nt,0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma, marker = '*', c='c',label =  'Log Mean')
    ax2[nt,1].scatter(df.log_abundance,df.log_sigma, c='c')
    ax2[nt,1].text(1.2,0.5,str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
    ax2[nt,0].semilogy()
    ax2[nt,1].semilogy()
    ax2[nt,1].semilogx()

    raw_r,raw_p  = scipy.stats.pearsonr(df['abundance'],df['sigma'])
    log_r,log_p = scipy.stats.pearsonr(df['log_abundance'],df['log_sigma'])
    print(raw_r,log_r)
    
    ax3[0,nt].hist(raw_r,color = 'red')
    ax3[0,0].text(0.3,0.5,'raw',horizontalalignment='center', verticalalignment='center', transform=ax3[nt,0].transAxes)
    ax3[1,nt].hist(log_r,color = 'b')
    ax3[1,0].text(0.3,0.5,'logged',horizontalalignment='center', verticalalignment='center', transform=ax3[nt,1].transAxes)


#working of graph stretching
    ax3[nt,1].set_xlim((stretch_rat[0])*df.abundance[0]),(stretch_rat[0])*df.abundance[-1])

    
    #plt.xlim(0, best_ratio)
#plt.ylim(0, best_ratio)
l1 = ax1[0,0].legend(loc = 'upper left')
l1.draw_frame(False)
l2 = ax2[1,0].legend(loc = 'upper left')
l2.draw_frame(False)

'''

'''#set y lim and xlim'
#yrange/xrange' pick wichever ratio is bigger logged vs unlogged'''


#raw and logged (std vs mean)
#loggrange , logx range, raw y range, raw xrange, log ratio and raw ratio, then find np.max(log or raw), then get it to do somthing prety
#convert small to largge ratio via factor and then use that factor to stretch axes lim of smaller to equal the biger
#get yranges to match and get y and x lim print out all numbers to the graph that we are making








plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.show()
fig1.savefig('../figures/raw_abiotic_dynamics.png')
fig2.savefig('../figures/logged_abiotic_dynamics.png')
fig3.savefig('../figures/coorelations_abiotic.png')



plt.show











print("I'm done with Abiotic Spike Assays bro! ")

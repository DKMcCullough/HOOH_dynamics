'''

name:   model_spiked_abiotic_batch_2.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model of Monoculture BCC assays to graph 0 H phyotplankton biomass and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop of all dfs in df_all for model and intits? 

'''


#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ODElib
import random as rd
import sys


######################################################
#reading in data and configureing 
#####################################################
df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_2-5-dataset', header = 1)


#df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_a = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  

df = df_a
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['log4'] = np.log(df['rep4'])
df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
df['avg2'] = df[['rep2', 'rep4']].mean(axis=1)
df['abundance'] = df[['rep1','rep2','rep3', 'rep4']].mean(axis=1)
df['std1'] = df[['rep1', 'rep3']].std(axis=1)
df['std2'] = df[['rep2', 'rep4']].std(axis=1)
df['sigma'] = df[['rep1','rep2','rep3', 'rep4']].std(axis=1)

df['lavg1'] = df[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
df['lavg2'] = df[['log2', 'log4']].mean(axis=1)
df['log_abundance'] = df[['log1','log2', 'log3','log4']].mean(axis=1)
df['stdlog1'] = df[['log1', 'log3']].std(axis=1) #taking stdv of logged reps
df['stdlog2'] = df[['log2', 'log4']].std(axis=1)
df['log_sigma'] = df[['log1','log2', 'log3','log4']].std(axis=1)

df['log_sigma'] = 0.2
df.loc[df['organism'] == 'H', 'log_sigma'] = 0.08

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False)]  #assay 0 H 
df4 = df.loc[(df['assay'].str.contains('4', case=False))]


## Reading in inits files for 0 and 400 models respectively

inits4 = pd.read_csv("../data/inits/abiotic_spike_2.csv")
#priors4 = inits4.to_dict()

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################

#actual model that will be run by the odelib model framework
def abiotic(y,t,params):
    deltah,Sh = params[0], params[1]
    H = y[0]
    dHdt = Sh - deltah*H 
    return [dHdt]



#initiating the model as a class in odelib (give us use of the methods in this class - like integrate :-) 
def get_model(df,priors):
    a1=ODElib.ModelFramework(ODE=abiotic,
                          parameter_names=['deltah','Sh', 'H0'],
                          state_names = snames,
                          dataframe=df,
                          deltah = deltah_prior.copy(),
                          Sh = Sh_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          H = H0_mean
                         )
    return a1


#####################################################
#####################################################
#####################################################



#find closest time 
def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)


#####################################################
#model param and state variable set up 
#####################################################

# state variable names
snames = ['H']

#sigma we give model to search withi for each param
pw = 1

#setting param prior guesses and inititaing as an odelib param class in odelib
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':3})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':100})

priors = {'deltah' :deltah_prior,'Sh' : Sh_prior,'H0' : H0_prior} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = inits4['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 100000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a4 = get_model(df4,priors) 


# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod4 = a4.integrate()

a4res = get_residuals(a4)

#########################################################
# graphing df and models together
#########################################################
c0 = 'blueviolet'


# Set up graph for Dynamics and param histograms

fig1,ax1 = plt.subplots(1,3,figsize=[10,7]) #plot creation and config 
#set titles of subplots
fig1.suptitle('Abiotic HOOH Model Output',fontsize = '16') #full title config
fig1.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30) #shift white space for better fig view
ax1[0].set_title('HOOH Dynamics',fontsize = '12')
ax1[0].set_ylabel('HOOH Concentration nM/mL',fontsize = '12')
ax1[0].set_xlabel('Time (days)',fontsize = '12')
ax1[1].set_title('Sh',fontsize = '12')
ax1[1].set_ylabel('Frequency',fontsize = '12')
ax1[1].set_xlabel('Parameter Value',fontsize = '12')
ax1[2].set_title('deltah',fontsize = '12')
ax1[2].set_ylabel('Frequency',fontsize = '12')
ax1[2].set_xlabel('Parameter Value',fontsize = '12')



#plot dynamics of data and model for 0 assay 
ax1[0].plot(df4.time,df4.abundance, marker='o',color = c0, label = 'abiotic - 4 H ') #data of 0 H assay
ax1[0].errorbar(df4.time,df4.abundance, yerr = df4.sigma, marker='o',color = c0) #data of 0 H assay
ax1[0].plot(mod4.time,mod4['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a4.plot_uncertainty(ax1[0],posteriors4,'H',100)

# plot histograms of params next to dynamics graphs
ax1[1].hist((posteriors4.Sh), facecolor=c0) #graphing Sh of 0 H assay 
ax1[2].hist((posteriors4.deltah), facecolor=c0) #graphing deltah of 0 H assay 

#config legends
l1 = ax1[0].legend(loc = 'lower right')
l1.draw_frame(False)




fig1.savefig('../figures/abiotic_4_dynamics')

########################################
#graph parameters against one another 
########################################

#graph set up

fig2,ax2 = plt.subplots(1,2, figsize=[9,6])
fig2.suptitle('Parameter Interactions ',fontsize = '14')

ax2[0].set_ylabel('deltah',fontsize = '12')
ax2[0].set_xlabel('Sh',fontsize = '12')
ax2[1].set_ylabel('log (deltah)',fontsize = '12')
ax2[1].set_xlabel('log (Sh)',fontsize = '12')

#ax2.set_title('0 HOOH')

plt.legend()
#adding text for more labels of graph

fig2.subplots_adjust(right=0.90, left=0.15,wspace = 0.45, hspace = 0.30) #shift white space for better fig view

#graphing each assay's parameters against each other 
ax2[0].scatter(posteriors4.Sh,posteriors4.deltah,color = c0)
ax2[1].scatter(np.log(posteriors4.Sh),np.log(posteriors4.deltah),color = c0)

#ax2[1,0].set_yscale('log')


#show full graph and save fig

fig2.savefig('../figures/abiotic_4_params')


#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig3,ax3 = plt.subplots(1,2,sharex=True,figsize=[8,4]) #make plot
fig3.suptitle('Trace plots for Logged Params ',fontsize = '16') #set main title 
fig3.subplots_adjust(right=0.90, wspace = 0.45, top = 0.85) #shift white space for better fig view
fig3.supxlabel('Model Iteration', fontsize = '14') #set overall x title 
#ax3[0].set_title('0 HOOH')
#ax3[1].set_title('400 HOOH ')
ax3[0].set_ylabel('Sh',fontsize = '12')
ax3[1].set_ylabel('deltah',fontsize = '12')

#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax3[0].scatter(posteriors4.iteration,(posteriors4.Sh),color = c0)
ax3[1].scatter(posteriors4.iteration,(posteriors4.deltah),color = c0)



#print out plot
fig3.savefig('../figures/abiotic_4_TRACE')


#########################################
#graphing Residuals of best model vs data 
##########################################

#making and confing of residuals plot
fig4, (ax0,ax1)= plt.subplots(1,2,figsize = (8,4)) #fig creationg of 1 by 2
fig4.suptitle('Abiotic HOOH 400 spike Model',fontsize = '16') #setting main title of fig

####### fig config and naming 

fig4.subplots_adjust(right=0.9, wspace = 0.45, hspace = 0.20, top = 0.8)

ax0.semilogy()
ax0.set_title('HOOH dynamics ',fontsize = '16')
ax1.set_title('Model residuals',fontsize = '12')

ax0.set_xlabel('Time (days)',fontsize = '12')
ax0.set_ylabel('HOOH (nM)',fontsize = '12')
ax1.set_ylabel('Data H value',fontsize = '12')
ax1.set_xlabel('Residual',fontsize = '12')



#model and residuals
ax0.plot(df4.time,df4.abundance, marker='o',color = c0, label = 'HOOH data mean') #data of 0 H assay
ax0.errorbar(df4.time,df4.abundance, yerr = df4.sigma, marker='o',color = c0) #data of 0 H assay
ax0.plot(mod4.time,mod4['H'],c='r',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a4.plot_uncertainty(ax0,posteriors4,'H',100)

ax1.scatter(a4res['res'], a4res['abundance'],color = c0, label = '400H spike')

#printing off graph
l4 = ax0.legend(loc = 'upper right')
l4.draw_frame(False)

plt.show()


fig4.savefig('../figures/abiotic_4_residuals')


#saving best params as the new inits file for the next file run. 

pframe4 = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe4.to_csv("../data/inits/abiotic_spike_2.csv")


# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')






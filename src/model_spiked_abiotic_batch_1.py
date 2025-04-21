'''

name:   model_spiked_abiotic_batch_1.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model of Monoculture BCC assays to graph 0 H phyotplankton biomass and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop of all dfs in df_all for model and intits? 

'''


#read in needed packages 
import helpers as hp
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

# get data and visualize uncertainty
df = hp.get_data('abiotic')
df4 = df.loc[(df['assay'].str.contains('4', case=False))]
sigma4 = hp.get_uncertainty(df4)
df4.loc[df['organism'] == 'H', 'log_sigma'] = sigma4
figa,axa = hp.plot_uncertainty(df4,sigma4)

## Reading in inits files for 0 and 400 models respectively
inits4 = pd.read_csv("../data/inits/abiotic_spike_1.csv")

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
def get_model(df):
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
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':360})


#setting H mean for odelib search 
H0_mean = inits4['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 10000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a4 = get_model(df4) 

# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod4 = a4.integrate()
a4res = get_residuals(a4)

#########################################################
# graphing df and models together
#########################################################
c0 = 'slateblue'


# Set up graph for Dynamics and param histograms

fig1,ax1 = plt.subplots(1,3,figsize=[12,4]) #plot creation and config 
fig1.subplots_adjust(wspace=0.3,bottom=0.2) #shift white space for better fig view
ax1[0].set_ylabel(r'H$_2$O$_2$ (Concentration pmol mL$^{-1}$)', fontsize = 12)
ax1[0].set_xlabel('Time (days)', fontsize = 12)
ax1[1].set_ylabel('Frequency', fontsize = 12)
ax1[1].set_xlabel('H$_2$O$_2$ supply rate \n ($S_H$, pmol mL$^{-1}$ day$^{-1}$)', fontsize = 12)
ax1[2].set_ylabel('Frequency', fontsize = 12)
ax1[2].set_xlabel('H$_2$O$_2$ decay rate \n ($\delta_H$, day$^{-1}$)', fontsize = 12)

#ax1[0].set_ylim([20, 600])

for (ax,l) in zip(ax1,'abc'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#plot dynamics of data and model for 0 assay 
ax1[0].plot(df4.time,df4.abundance, marker='o',color = c0, label = r'H$_2$O$_2$ data ') #data of 0 H assay
ax1[0].plot(mod4.time,mod4['H'],c='k',lw=1.5,label='Model best fit') #best model fit of 0 H assay
ax1[0].errorbar(df4.time,df4.abundance, yerr = df4.sigma, marker='o',color = c0) #data of 0 H assay
a4.plot_uncertainty(ax1[0],posteriors4,'H',100)

# plot histograms of params next to dynamics graphs
ax1[1].hist((posteriors4.Sh), facecolor=c0) #graphing Sh of 0 H assay 
ax1[2].hist((posteriors4.deltah), facecolor=c0) #graphing deltah of 0 H assay 

#config legends
l1 = ax1[0].legend(loc = 'lower right')
l1.draw_frame(False)


fig1.savefig('../figures/abiotic1_4_dynamics')

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

fig2.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view

#graphing each assay's parameters against each other 
ax2[0].scatter(posteriors4.Sh,posteriors4.deltah,color = c0)
ax2[1].scatter(np.log(posteriors4.Sh),np.log(posteriors4.deltah),color = c0)

#ax2[1,0].set_yscale('log')


#show full graph and save fig

#fig2.savefig('../figures/abiotic1_4_params')


#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig3,ax3 = plt.subplots(1,2,sharex=True,figsize=[8,4]) #make plot
fig3.suptitle('Trace plots for H Params ',fontsize = '16') #set main title 
fig3.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view
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
#fig3.savefig('../figures/abiotic1_4_TRACE')


#########################################
#graphing Residuals of best model vs data 
##########################################

#making and confing of residuals plot
fig4, (ax0,ax1)= plt.subplots(1,2,figsize = (8,4)) #fig creationg of 1 by 2
fig4.suptitle('Abiotic HOOH 400 spike Model',fontsize = '16') #setting main title of fig

####### fig config and naming 

fig4.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2)

ax0.semilogy()
ax0.set_title('HOOH dynamics ',fontsize = '16')
ax1.set_title('Model residuals',fontsize = '12')

ax0.set_xlabel('Time (days)',fontsize = '12')
ax0.set_ylabel('HOOH (nM)',fontsize = '12')
ax1.set_ylabel('Data H value',fontsize = '12')
ax1.set_xlabel('Residual',fontsize = '12')

ax0.set_ylim([20, 600])
ax1.set_ylim([20, 600])


#model and residuals
ax0.plot(df4.time,df4.abundance, marker='o',color = c0, label = 'HOOH data mean') #data of 0 H assay
ax0.errorbar(df4.time,df4.abundance, yerr = df4.sigma, marker='o',color = c0) #data of 0 H assay
ax0.plot(mod4.time,mod4['H'],c='r',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a4.plot_uncertainty(ax0,posteriors4,'H',100)

ax1.errorbar(a4res['res'], a4res['abundance'],yerr=df4.sigma,color = c0,marker = 'o', markersize = 4, ls = 'none',elinewidth=2,label = '400H spike')

#printing off graph
l4 = ax0.legend(loc = 'lower right')
l4.draw_frame(False)

#fig4.savefig('../figures/abiotic1_4_residuals')


#saving best params as the new inits file for the next file run. 

pframe4 = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe4.to_csv("../data/inits/abiotic_spike_1.csv")


# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')






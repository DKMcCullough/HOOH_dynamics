'''

name:   model_abiotic_batch_1.py 

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
import helpers as hp
import ODElib
import random as rd
import sys

######################################################
#reading in data and configureing 
#####################################################
df = hp.get_data('abiotic')

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False)]  #assay 0 H 
sigma0 = hp.get_uncertainty(df0)
df0.loc[df['organism'] == 'H', 'log_sigma'] = sigma0
figa,axa = hp.plot_uncertainty(df0,sigma0)
#figa.savefig('../figures/error_0spike')

## Reading in inits files for 0 and 400 models respectively
inits0 = pd.read_csv("../data/inits/abiotic_control_1.csv")

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

#priors

#setting param prior guesses and inititaing as an odelib param class in odelib
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':20})

priors = {'deltah' :deltah_prior,'Sh' : Sh_prior,'H0' : H0_prior} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = inits0['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 100000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a0 = get_model(df0) 

# do fitting
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod0 = a0.integrate()

a0res = get_residuals(a0)  #is this using the best fit or just a first run???

#########################################################
# graphing df and models together
#########################################################

c0 = 'plum'


# Set up graph for Dynamics and param histograms

fig1,ax1 = plt.subplots(1,3,figsize=[12,4]) #plot creation and config 
#set titles of subplots
#fig1.suptitle('Abiotic HOOH Model Output', fontsize = 14) #full title config
fig1.subplots_adjust(wspace=0.3,bottom=0.2) #shift white space for better fig view
ax1[0].set_ylabel(r'H$_2$O$_2$ (Concentration pmol mL$^{-1}$)', fontsize = 12)
ax1[0].set_xlabel('Time (days)', fontsize = 12)
ax1[1].set_ylabel('Probability density', fontsize = 12)
ax1[1].set_xlabel('H$_2$O$_2$ supply rate \n ($S_H$, pmol mL$^{-1}$ day$^{-1}$)', fontsize = 12)
ax1[2].set_ylabel('Probability density', fontsize = 12)
ax1[2].set_xlabel('H$_2$O$_2$ decay rate \n ($\delta_H$, day$^{-1}$)', fontsize = 12)

#ax1[0].set_ylim([20, 600])

for (ax,l) in zip(ax1,'abc'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#plot dynamics of data and model for 0 assay 
ax1[0].plot(df0.time,df0.abundance, marker='o',color = c0, label = r'H$_2$O$_2$ data ') #data of 0 H assay
ax1[0].plot(mod0.time,mod0['H'],c='darkred',lw=1.5,label='Model best fit') #best model fit of 0 H assay

a0.plot_uncertainty(ax1[0],posteriors0,'H',100)

# plot histograms of params next to dynamics graphs
ax1[1].hist(((posteriors0.Sh)),density=True, facecolor=c0) #graphing Sh of 0 H assay 
ax1[2].hist(((posteriors0.deltah)),density=True, facecolor=c0) #graphing deltah of 0 H assay 

#config legends
l1 = ax1[0].legend(loc = 'lower right')
l1.draw_frame(False)

#save fig
fig1.savefig('../figures/abiotic1_0_dynamics',bbox_inches='tight')
fig1.savefig('../figures/figure2.tiff',bbox_inches='tight',dpi=300,format='tiff')

########################################
#graph parameters against one another 
########################################

#graph set up

fig2,ax2 = plt.subplots(1,2, figsize=[7,4])
fig2.suptitle('Parameter Interactions ')

ax2[0].set_ylabel('deltah', fontsize = 12)
ax2[0].set_xlabel('Sh', fontsize = 12)
ax2[1].set_ylabel('log (deltah)', fontsize = 12)
ax2[1].set_xlabel('log (Sh)', fontsize = 12)


#adding text for more labels of graph

fig2.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view

#graphing each assay's parameters against each other 
ax2[0].scatter(posteriors0.Sh,posteriors0.deltah,color = c0)
ax2[1].scatter(np.log(posteriors0.Sh),np.log(posteriors0.deltah),color = c0)
plt.legend()
#save fig
#fig2.savefig('../figures/abiotic1_0_params')


#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig3,ax3 = plt.subplots(1,2,sharex=True,figsize=[8,4]) #make plot
fig3.suptitle('Trace plots for H Params ', fontsize = 14) #set main title 
fig3.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view
fig3.supxlabel('Model Iteration', fontsize = 12) #set overall x title 

ax3[0].set_ylabel('Log Sh', fontsize = 12)
ax3[0].set_xlabel('Model iteration', fontsize = 12)
ax3[1].set_ylabel('Log deltah', fontsize = 12)
ax3[1].set_xlabel('Model iteration', fontsize = 12)
#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax3[0].scatter(posteriors0.iteration,posteriors0.Sh,color = c0)
ax3[1].scatter(posteriors0.iteration,posteriors0.deltah,color = c0)


#print out plot
#fig3.savefig('../figures/abiotic1_0_TRACE')

pframe0 = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe0.to_csv("../data/inits/abiotic_control_1.csv")

# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



'''

name:   model_spiked_abiotic_batch.py 

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
df_1 = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_1-31-dataset', header = 1)
df_2 = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_2-5-dataset', header = 1)

df_all = df_1


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

inits4 = pd.read_csv("../data/inits/abiotic4.csv")
#priors4 = inits4.to_dict()

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
'''
def set_best_params(model,posteriors,snames):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    bestchain = posteriors.iloc[im]["chain#"]
    posteriors = posteriors[posteriors["chain#"]==bestchain]
    model.set_parameters(**posteriors.loc[im][a4.get_pnames()].to_dict())
    model.set_inits(**{o:posteriors.loc[im][a4.get_pnames()].to_dict()[o+'0'] for o in ['H']})

###############
#####only set for 0 a for idk if 400 model is working correctly. #######

#function for plotting uncertainty once model has been run 
def plot_uncertainty(ax,model,posteriors,ntimes):
    for a in range(ntimes):
        im = rd.choice(posteriors.index) 
        model.set_inits(**{'H':posteriors.loc[im][model.get_pnames()].to_dict()['H0']})
        model.set_parameters(**posteriors.loc[im][model.get_pnames()].to_dict())
        mod = model.integrate()
        ax.plot(mod.time,mod['H'],c=str(0.8),lw=1,zorder=1)

'''
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
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':300})

priors = {'deltah' :deltah_prior,'Sh' : Sh_prior,'H0' : H0_prior} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = inits4['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 10000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a4 = get_model(df4,priors) 

#broken here!!!!!!!!!!
# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()

a4res = get_residuals(a4)

#########################################################
# graphing df and models together
#########################################################
c0 = 'darkviolet'


# Set up graph for Dynamics and param histograms

fig1,ax1 = plt.subplots(1,3,figsize=[10,7]) #plot creation and config 
#set titles of subplots
fig1.suptitle('Abiotic HOOH Model Output') #full title config
fig1.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30) #shift white space for better fig view
ax1[0].set_title('400 H Model-Data Dynamics')
fig1.supylabel('HOOH Concentration')
ax1[1].set_title('Sh')
ax1[2].set_title('deltah')
ax1[1].text(2.3, -130, 'Frequency')


fig1.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30) #shift white space for better fig view

#config legends
l1 = ax1[0].legend(loc = 'lower right')
l2 = ax1[1].legend(loc = 'upper left')
l1.draw_frame(False)
l2.draw_frame(False)

#graph dynamics of data and model (best model line), and posterior guesses (grey lines) of 0 and 400 respectively

#plot dynamics of data and model for 0 assay 
ax1[0].plot(df4.time,df4.abundance, marker='o',color = c0, label = 'abiotic - 4 H ') #data of 0 H assay
ax1[0].plot(mod4.time,mod4['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a4.plot_uncertainty(ax1[0],posteriors4,'H',100)
#plot_uncertainty(ax1[0],a4,posteriors4,100) #plotting 100 itterations of model search for 0 H assay 
#plot 400 assay dynamics and models
#ax1[1].plot(df4.time,df4.abundance, marker='o',color = c4, label = 'abiotic - 400 H ')#data of 400 H

# plot histograms of params next to dynamics graphs
ax1[1].hist((np.log(posteriors4.Sh)), facecolor=c0) #graphing Sh of 0 H assay 
ax1[2].hist((np.log(posteriors4.deltah)), facecolor=c0) #graphing deltah of 0 H assay 

fig1.savefig('../figures/abiotic_4_dynamics')

########################################
#graph parameters against one another 
########################################

#graph set up

fig2,ax2 = plt.subplots(1,2, figsize=[9,6])
fig2.suptitle('Parameter Interactions ')

ax2[0].set_ylabel('deltah')
ax2[0].set_xlabel('Sh')
ax2[1].set_ylabel('ln (deltah)')
ax2[1].set_xlabel('ln (Sh)')

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
fig3,ax3 = plt.subplots(1,2,sharex=True,figsize=[10,5]) #make plot
fig3.suptitle('Trace plots for Logged Params ') #set main title 
fig3.subplots_adjust(right=0.90, wspace = 0.45, top = 0.85) #shift white space for better fig view
fig3.supxlabel('Model Iteration') #set overall x title 
#ax3[0].set_title('0 HOOH')
#ax3[1].set_title('400 HOOH ')
ax3[0].set_ylabel('Log Sh')
ax3[1].set_ylabel('Log deltah')

#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax3[0].scatter(posteriors4.iteration,np.log(posteriors4.Sh),color = c0)
ax3[1].scatter(posteriors4.iteration,np.log(posteriors4.deltah),color = c0)



#print out plot
fig3.savefig('../figures/abiotic_4_TRACE')


#########################################
#graphing Residuals of best model vs data 
##########################################

#making and confing of residuals plot
fig4,ax4 = plt.subplots(figsize=[8,5])
fig4.suptitle('Residuals vs Model Fit Value ')
fig4.supylabel('Model Value (H)')
fig4.supxlabel('Residual')
#config legends for data differentialtion 
l4 = ax4.legend()
l4.draw_frame(False)

#plotting residual function output residual and abundance columns 
ax4.scatter(a4res['res'], a4res['abundance'],label = '0 H', color = c0) #where )

#how to get residuals from all posterior runs not just best???

#print out plot
fig4.savefig('../figures/abiotic_4_residuals')

pframe4 = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe4.to_csv("../data/inits/abiotic4.csv")


# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print('\n Done my guy \n')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')






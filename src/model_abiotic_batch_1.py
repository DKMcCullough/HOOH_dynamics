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
import ODElib
import random as rd
import sys

plt.rcParams["font.family"] = "Times New Roman"

######################################################
#reading in data and configureing 
#####################################################
df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_1-31-dataset', header = 1)




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
nits = 10000


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

fig1,ax1 = plt.subplots(1,3,figsize=[12,5]) #plot creation and config 
#set titles of subplots
#fig1.suptitle('Abiotic HOOH Model Output', fontsize = 14) #full title config
fig1.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view
ax1[0].set_title(r'H$_2$O$_2$ Dynamics', fontsize = 12)
ax1[0].set_ylabel(r'H$_2$O$_2$ Concentration nM/mL', fontsize = 12)
ax1[0].set_xlabel('Time (days)', fontsize = 12)
ax1[1].set_title(r'$S_H$', fontsize = 12)
ax1[1].set_ylabel('Frequency', fontsize = 12)
ax1[1].set_xlabel('Parameter Value', fontsize = 12)
ax1[2].set_title(r'$\delta_H$', fontsize = 12)
ax1[2].set_ylabel('Frequency', fontsize = 12)
ax1[2].set_xlabel('Parameter Value', fontsize = 12)

#ax1[0].set_ylim([20, 600])

for (ax,l) in zip(ax1,'abc'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#plot dynamics of data and model for 0 assay 
ax1[0].plot(df0.time,df0.abundance, marker='o',color = c0, label = r'H$_2$O$_2$ data ') #data of 0 H assay
ax1[0].plot(mod0.time,mod0['H'],c='darkred',lw=1.5,label='Model best fit') #best model fit of 0 H assay

a0.plot_uncertainty(ax1[0],posteriors0,'H',100)

# plot histograms of params next to dynamics graphs
ax1[1].hist(((posteriors0.Sh)), facecolor=c0) #graphing Sh of 0 H assay 
ax1[2].hist(((posteriors0.deltah)), facecolor=c0) #graphing deltah of 0 H assay 

#config legends
l1 = ax1[0].legend(loc = 'upper right')
l1.draw_frame(False)

#save fig
fig1.savefig('../figures/abiotic1_0_dynamics')

plt.show()

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
fig2.savefig('../figures/abiotic1_0_params')


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
fig3.savefig('../figures/abiotic1_0_TRACE')


#########################################
#graphing Residuals of best model vs data 
##########################################

fig4, (ax0,ax1)= plt.subplots(1,2,figsize = (8,4)) #fig creationg of 1 by 2
fig4.suptitle('Abiotic HOOH Model',fontsize = '14') #setting main title of fig

fig4.subplots_adjust(right=0.9, wspace = 0.45, hspace = 0.20)

ax0.semilogy()
ax0.set_title('HOOH dynamics ',fontsize = '14')
ax1.set_title('Model residuals',fontsize = '14')

ax0.set_xlabel('Time (days)', fontsize = 12)
ax0.set_ylabel('HOOH (nM)', fontsize = 12)
ax1.set_ylabel('HOOH data value', fontsize = 12)
ax1.set_xlabel('Residual', fontsize = 12)

ax0.set_ylim([20, 600])
ax1.set_ylim([20, 600])



#model and residuals
ax0.plot(df0.time,df0.abundance, marker='o',color = c0, label = 'HOOH data mean') #data of 0 H assay
ax0.errorbar(df0.time,df0.abundance, yerr = df0.sigma, marker='o',color = c0) #data of 0 H assay
ax0.plot(mod0.time,mod0['H'],c='darkred',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a0.plot_uncertainty(ax0,posteriors0,'H',100)

ax1.errorbar(a0res['res'], a0res['abundance'],yerr=df4.sigma,color = c0,marker = 'o', markersize = 4, ls = 'none',elinewidth=2,label = '0H spike')

#printing off graph
l4 = ax0.legend(loc = 'lower right')
l4.draw_frame(False)

plt.show()
fig4.savefig('../figures/abiotic1_0_residuals')




pframe0 = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe0.to_csv("../data/inits/abiotic_control_1.csv")



# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



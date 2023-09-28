'''

name:   model_abiotic_batch.py 

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


#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
#plt.rcParams.update({'font.size': 12})
#plt.rcParams['legend.fontsize'] = 'small'


#####################################################
# read in data and formatting
#####################################################

#main df read in 
df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
#df_all = pd.read_csv("../data/BCC_2-5-dataset.csv",header=1)

#format empty columns and column names 
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(hr)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#split df into only abiotice H data 
df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  

# configuring df for modeling via odelib
df = df_abiotic  #ensuring df is working df
df['log_abundance'] = np.log(df['raw_abundance']) #creating logabundance value and column from raw df data column 
#df['log_sigma'] = np.std(df['HOOH_stdv']) #maybe use eventually once we have actual replicates 
df['log_sigma'] = 0.04 # made up number to have a stdv for our model to fit 
df = df.rename(columns={"raw_abundance": "abundance"}) #renaming raw to abundance for odelib to graph against

#splitting df into 0 HOOH and 400 HOOH assay dfs 
df0 = df.loc[~ df['assay'].str.contains('4', case=False)] 
df4 = df.loc[df['assay'].str.contains('4', case=False)] 

## Reading in inits files for 0 and 400 models respectively
inits0 = pd.read_csv("../data/inits/abiotic0.csv")
inits4 = pd.read_csv("../data/inits/abiotic4.csv")
#inits0 = pd.read_csv(("../data/inits/"+str(s)+"_abiotic0".csv"))
#inits = pd.read_csv(("../data/inits/"+str(s)+"_abiotic"+str(t)+".csv"))
t = [0,4]

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
                          H = H0_mean,
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
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})

Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':4})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':360})

#setting H mean for odelib search 
H0_mean = df.loc[df['time'] == 0, 'abundance'].iloc[0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 1000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a0 = get_model(df0) # initialising model for 0 df
a4 = get_model(df4) # initialising model for 400 df
 

# do fitting for 0 an 400 model 
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1, print_report=True)
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1, print_report=True)


# set best params for 0 an 400 model 
#set_best_params(a0,posteriors0,snames)
#set_best_params(a4,posteriors4,snames)
a0.set_best_params(posteriors0)   #new func working for best params 
a4.set_best_params(posteriors4)   #new func workinng  for best params 

# run model with optimal params for 0 an 400 model 
mod0 = a0.integrate()
mod4 = a4.integrate()

#get residuals from model 
a0res = get_residuals(a0)  #is this using the best fit or just a first run???
a4res = get_residuals(a4)


#########################################################
# graphing df and models together
#########################################################
c0 = 'darkgreen'
c4 = 'darkviolet'

# Set up graph for Dynamics and param histograms

fig1,ax1 = plt.subplots(2,3,figsize=[10,7]) #plot creation and config 
#set titles of subplots
fig1.suptitle('Abiotic HOOH Model Output') #full title config
fig1.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30) #shift white space for better fig view
ax1[0,0].set_title('Model-Data Dynamics')
fig1.supylabel('HOOH Concentration nM/mL')
ax1[0,1].set_title('Sh')
ax1[0,2].set_title('deltah')
ax1[1,1].text(3.5, -15, 'Frequency')


fig1.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30) #shift white space for better fig view

#config legends
#l1 = ax1[0,0].legend(loc = 'lower right')
#l2 = ax1[1,0].legend(loc = 'upper left')
#l1.draw_frame(False)
#l2.draw_frame(False)

#graph dynamics of data and model (best model line), and posterior guesses (grey lines) of 0 and 400 respectively

#plot dynamics of data and model for 0 assay 
ax1[0,0].plot(df0.time,df0.abundance, marker='o',color = c0, label = 'abiotic - 0 H ') #data of 0 H assay
ax1[0,0].plot(mod0.time,mod0['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
#plot_uncertainty(ax1[0,0],a0,posteriors0,100) #plotting 100 itterations of model search for 0 H assay 
a0.plot_uncertainty(ax1[0,0],posteriors0,'H',100) #new uncertainty function working 
#plot 400 assay dynamics and models
ax1[1,0].plot(df4.time,df4.abundance, marker='o',color = c4, label = 'abiotic - 400 H ')#data of 400 H
ax1[1,0].plot(mod4.time,mod4['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 400 H assay 
#plot_uncertainty(ax1[1,0],a4,posteriors4,100) #plotting 100 itterations of model for 400 H assay 
a4.plot_uncertainty(ax1[1,0],posteriors4,'H',100)  #new uncertainty function working 
# plot histograms of params next to dynamics graphs
ax1[0,1].hist((np.log(posteriors0.Sh)), facecolor=c0) #graphing Sh of 0 H assay 
ax1[0,2].hist((np.log(posteriors0.deltah)), facecolor=c0) #graphing deltah of 0 H assay 
ax1[1,1].hist((np.log(posteriors4.Sh)),facecolor=c4) #graphing Sh of 400 H assay 
ax1[1,2].hist((np.log(posteriors4.deltah)),facecolor=c4) #graphing deltah of 400 H assay 

fig1.savefig('../figures/abiotic_0and400_dynamics')

########################################
#graph parameters against one another 
########################################

#graph set up

fig2,ax2 = plt.subplots(2,2, figsize=[8,5])
fig2.suptitle('Deltah vs Sh ')
fig2.supylabel('ln deltah')
fig2.supxlabel('ln Sh')
ax2[0,0].set_title('0 HOOH')
ax2[0,1].set_title('400 HOOH ')
#plt.legend()
#adding text for more labels of graph
ax2[0,1].text(22.5, 0.01, 'RAW',)
ax2[1,1].text(22.5, (-4.8), 'LOG',)
fig2.subplots_adjust(right=0.90, left=0.15,wspace = 0.25, hspace = 0.30) #shift white space for better fig view

#graphing each assay's parameters against each other 
ax2[0,0].scatter(posteriors0.Sh,posteriors0.deltah,color = c0)
ax2[0,1].scatter(posteriors4.Sh,posteriors4.deltah,color = c4)
ax2[1,0].scatter(np.log(posteriors0.Sh),np.log(posteriors0.deltah),color = c0)
ax2[1,1].scatter(np.log(posteriors4.Sh),np.log(posteriors4.deltah),color = c4)

#ax2[1,0].set_yscale('log')


#show full graph and save fig

fig2.savefig('../figures/abiotic_0and400_params')


#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig3,ax3 = plt.subplots(2,2,sharex=True,figsize=[8,5]) #make plot
fig3.suptitle('Trace plots for Logged Params ') #set main title 
fig3.subplots_adjust(right=0.90, wspace = 0.25, top = 0.85) #shift white space for better fig view
fig3.supxlabel('Model Iteration') #set overall x title 
ax3[0,0].set_title('0 HOOH')
ax3[0,1].set_title('400 HOOH ')
ax3[0,0].set_ylabel('Log Sh')
ax3[1,0].set_ylabel('Log deltah')

#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax3[0,0].scatter(posteriors0.iteration,np.log(posteriors0.Sh),color = c0)
ax3[0,1].scatter(posteriors4.iteration,np.log(posteriors4.Sh),color = c4)
ax3[1,0].scatter(posteriors0.iteration,np.log(posteriors0.deltah),color = c0)
ax3[1,1].scatter(posteriors4.iteration,np.log(posteriors4.deltah),color = c4)



#print out plot
fig3.savefig('../figures/abiotic_0and400_TRACE')


#########################################
#graphing Residuals of best model vs data 
##########################################

#making and confing of residuals plot
fig4,ax4 = plt.subplots(2,1,sharex = True,figsize=[8,5])
fig4.suptitle('Residuals vs Model Fit Value ')
fig4.supylabel('Model Value (H)')
fig4.supxlabel('Residual')
#config legends for data differentialtion 
#l4 = ax4[0].legend()
#l5 = ax4[1].legend()
#l4.draw_frame(False)
#l5.draw_frame(False)

#plotting residual function output residual and abundance columns 
ax4[0].scatter(a0res['res'], a0res['abundance'],label = '0 H', color = c0) #where )
ax4[1].scatter(a4res['res'], a4res['abundance'],label = '400 H', color = c4)

#how to get residuals from all posterior runs not just best???

#print out plot
fig4.savefig('../figures/abiotic_0and400_residuals')

plt.show()
# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print('\n Done my guy \n')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



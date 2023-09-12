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
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'


#####################################################
# read in data and formatting
#####################################################

#main df read in 
df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
#format empty columns and column names 
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#split df into only abiotice H data 
df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  


# configuring df for modeling via odelib
df = df_abiotic  #ensuring df is working df
df['log_abundance'] = np.log(df['raw_abundance']) #creating logabundance value and column from raw df data column 
#df['log_sigma'] = np.std(df['HOOH_stdv']) #maybe use eventually once we have actual replicates 
df['log_sigma'] = 0.1 # made up number to have a stdv for our model to fit 
df = df.rename(columns={"raw_abundance": "abundance"}) #renaming raw to abundance for odelib to graph against

#splitting df into 0 HOOH and 400 HOOH assay dfs 
df0 = df.loc[~ df['assay'].str.contains('4', case=False)] 
df4 = df.loc[df['assay'].str.contains('4', case=False)] 

## Reading in inits files for 0 and 400 models respectively
inits0 = pd.read_csv("../data/inits/abiotic0.csv")
inits4 = pd.read_csv("../data/inits/abiotic4.csv")


#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def set_best_params(model,posteriors,snames):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    bestchain = posteriors.iloc[im]["chain#"]
    posteriors = posteriors[posteriors["chain#"]==bestchain]
    model.set_parameters(**posteriors.loc[im][a0.get_pnames()].to_dict())
    model.set_inits(**{o:posteriors.loc[im][a0.get_pnames()].to_dict()[o+'0'] for o in ['H']})
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


#attatch time  
#find closesst time 
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

Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':6})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':(pw/10),'scale':1e+5})

#setting H mean for odelib search 
H0_mean = df.loc[df['time'] == 0, 'abundance'].iloc[0]

# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 100000

#has P, N H, if not OB - trie 

# get_models for 0 an 400 model 
a0 = get_model(df0) 
a4 = get_model(df4) 
 

# do fitting for 0 an 400 model 
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1, print_report=False)
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1, print_report=False)


# set best params for 0 an 400 model 
set_best_params(a0,posteriors0,snames)
set_best_params(a4,posteriors4,snames)

# run model with optimal params for 0 an 400 model 
mod0 = a0.integrate()
mod4 = a4.integrate()

#get residuals from model 
a0res = get_residuals(a0)  #is this using the best fit or just a first run???
a4res = get_residuals(a4)


#########################################################
# graphing df and models together
#########################################################

# Set up graph for Dynamics and param histograms

fig4,ax4 = plt.subplots(2,3,figsize=[10,7]) #plot creation and config 
fig4.suptitle('Abiotic HOOH Model Output') #full title config
fig4.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30) #shift white space for better fig view

#graph dynamics of data and model (best model line), and posterior guesses (grey lines) of 0 and 400 respectively

#plot dynamics of data and model for 0 and 400 assays  
ax4[0,0].plot(df0.time,df0.abundance, marker='o',label = 'abiotic - 0 H ') #data of 0 H assay
ax4[0,0].plot(mod0.time,mod0['H'],c='r',lw=1.5,label=' model best fit') #best model fit of 0 H assay
plot_uncertainty(ax4[0,0],a0,posteriors0,100) #plotting 100 itterations of model search for 0 H assay 
ax4[1,0].plot(df4.time,df4.abundance, marker='o',label = 'abiotic - 400 H ')#data of 400 H
ax4[1,0].plot(mod4.time,mod4['H'],c='r',lw=1.5,label=' model best fit') #best model fit of 400 H assay 
plot_uncertainty(ax4[1,0],a4,posteriors4,100) #plotting 100 itterations of model for 400 H assay 


# plot histograms of params next to dynamics graphs
ax4[0,1].hist((np.log(posteriors0.Sh))) #graphing Sh of 0 H assay 
ax4[0,2].hist((np.log(posteriors0.deltah))) #graphing deltah of 0 H assay 
ax4[1,1].hist((np.log(posteriors4.Sh))) #graphing Sh of 400 H assay 
ax4[1,2].hist((np.log(posteriors4.deltah))) #graphing deltah of 400 H assay 

fig4.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30) #shift white space for better fig view

#set titles of subplots
ax4[0,0].set_title('Model-Data Dynamics')
ax4[0,1].set_title('Sh')
ax4[0,2].set_title('deltah')
ax4[0,0].set_ylabel('HOOH Concentration nM/mL')
ax4[1,0].set_ylabel('HOOH Concentration nM/mL')

#config legend 
l1 = ax4[0,0].legend(loc = 'upper left')
l2 = ax4[1,0].legend(loc = 'upper left')
l1.draw_frame(False)
l2.draw_frame(False)

plt.show()

'''
#fig4.savefig('../figures/abiotic_odelib')

########################################
#graph parameters against one another 
########################################

#graph set up

fig5,ax5 = plt.subplots(2,1,sharex=True, figsize=[8,5])
fig5.suptitle('deltah vs Sh ')
#graphing each assay's parameters against each other 
ax5[0].scatter(posteriors0.Sh,posteriors0.deltah)
ax5[1].scatter(posteriors4.Sh,posteriors4.deltah)
ax5[0].set_yscale('log')
ax5[1].set_xscale('log')
#ax5[0].set_xlabel('Frequency Sh')
ax5[1].set_xlabel('Frequency Sh')
ax5[0].set_ylabel('Frequency deltah')
ax5[1].set_ylabel('Frequency deltah')
#ax5[0].set_yscale('log')
plt.legend()
plt.show()

#fig5.savefig('../figures/abiotic_params')
'''

#################################
#graphing logged parameter values
##################################
#crating and config of fig 6
fig6,ax6 = plt.subplots(2,2,sharex=True,figsize=[8,5]) #make plot
fig6.suptitle('Trace plots for Logged deltah and Sh ') #set main title 
fig6.subplots_adjust(right=0.90, wspace = 0.25, top = 0.85) #shift white space for better fig view
fig6.supxlabel('Model Iteration') #set overall x title 
ax6[0,0].set_title('0 HOOH')
ax6[0,1].set_title('400 HOOH ')
ax6[0,0].set_ylabel('log Sh')
ax6[1,0].set_ylabel('log deltah')

#graphing iteration number vs parameter numbert logged 
ax6[0,0].scatter(posteriors0.iteration,np.log(posteriors0.Sh))
ax6[0,1].scatter(posteriors4.iteration,np.log(posteriors4.Sh))
ax6[1,0].scatter(posteriors0.iteration,np.log(posteriors0.deltah))
ax6[1,1].scatter(posteriors4.iteration,np.log(posteriors4.deltah))



#print out plot
plt.show()

#########################################
#graphing Residuals of best model vs data 
##########################################

#making and confing of residuals plot
fig7,ax7 = plt.subplots(2,1,sharex = True,figsize=[8,5])
fig7.suptitle('Residuals vs Fit Value ')
fig7.supylabel('Model Value (H)')

fig7.supxlabel('Residual')

#plotting residual function output residual and abundance columns 
ax7[0].scatter(a0res['res'], a0res['abundance'],label = '0 H') #where )
ax7[1].scatter( a4res['res'],a4res['abundance'],label = '400 H')


#config legends for data differentialtion 
l7 = ax7[0].legend()
l8 = ax7[1].legend()
l7.draw_frame(False)
l8.draw_frame(False)

#print out plot
plt.show()


# 'program finished' flag
print('\n Done my guy \n')

print('\n Im free Im free! Im done calculating!' )




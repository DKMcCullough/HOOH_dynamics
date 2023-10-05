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


#####################################################
# read in data and formatting
#####################################################

df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df = df_all

#make unlogged avgs
df['avg1'] = df[['rep1', 'rep3']].mean(axis=1)
df['std1'] = df[['rep1', 'rep3']].std(axis=1)

#make logs of data
df['log1'] = np.log(df['rep1'])
df['log3'] = np.log(df['rep3'])

#avg and std of logged data
df['lavg1'] = df[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
df['stdlog1'] = df[['log1', 'log3']].std(axis=1) #taking stdv of logged reps

#splicing abiotic and mono or coculture data
df_abiotic = df.loc[df['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df.loc[df['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df.loc[~df['assay'].str.contains('coculture', case=False)].copy()  
df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 

#setting working df
df = df_abiotic

df = df.rename(columns={"avg1": "abundance"}) #renaming raw to abundance for odelib to graph against
df = df.rename(columns={"stdlog1": "log_sigma"}) #renaming raw to abundance for odelib to graph against


#splitting df into 0 HOOH and 400 HOOH assay dfs 
df0 = df.loc[~ df['assay'].str.contains('4', case=False)] 
df4 = df.loc[df['assay'].str.contains('4', case=False)] 

## Reading in inits files for 0 and 400 models respectively
inits0 = pd.read_csv("../data/inits/abiotic0.csv")
#priors0 = inits0.to_dict()

inits4 = pd.read_csv("../data/inits/abiotic4.csv")
#priors4 = inits4.to_dict()

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
def get_model(df,priors):
    a1=ODElib.ModelFramework(ODE=abiotic,
                          parameter_names=['deltah','Sh', 'H0'],
                          state_names = snames,
                          dataframe=df,
                          deltah = priors['deltah'],
                          Sh = priors['Sh'],
                          H0  = priors['H0'],
                          t_steps=1000,
                          H = priors['H0'],
                         )
    return a1
#vroken rn  - get priors frm 
 
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
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})

priors = {} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = df.loc[df['time'] == 0, 'abundance'].iloc[0]



# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 1000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a0 = get_model(df0,priors) # initialising model for 0 df
a4 = get_model(df4,priors) # initialising model for 400 df
 

# do fitting for 0 an 400 model 
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1, print_report=True)
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1, print_report=False)

#new way to use set best params funct
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
l1 = ax1[0,0].legend(loc = 'lower right')
l2 = ax1[1,0].legend(loc = 'upper left')
l1.draw_frame(False)
l2.draw_frame(False)

#graph dynamics of data and model (best model line), and posterior guesses (grey lines) of 0 and 400 respectively

#plot dynamics of data and model for 0 assay 
ax1[0,0].plot(df0.time,df0.abundance, marker='o',color = c0, label = 'abiotic - 0 H ') #data of 0 H assay
ax1[0,0].plot(mod0.time,mod0['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
plot_uncertainty(ax1[0,0],a0,posteriors0,100) #plotting 100 itterations of model search for 0 H assay 
#plot 400 assay dynamics and models
ax1[1,0].plot(df4.time,df4.abundance, marker='o',color = c4, label = 'abiotic - 400 H ')#data of 400 H
ax1[1,0].plot(mod4.time,mod4['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 400 H assay 
plot_uncertainty(ax1[1,0],a4,posteriors4,100) #plotting 100 itterations of model for 400 H assay 
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
plt.legend()
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
l4 = ax4[0].legend()
l5 = ax4[1].legend()
l4.draw_frame(False)
l5.draw_frame(False)

#plotting residual function output residual and abundance columns 
ax4[0].scatter(a0res['res'], a0res['abundance'],label = '0 H', color = c0) #where )
ax4[1].scatter(a4res['res'], a4res['abundance'],label = '400 H', color = c4)

#how to get residuals from all posterior runs not just best???

#print out plot
fig4.savefig('../figures/abiotic_0and400_residuals')


# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print('\n Done my guy \n')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



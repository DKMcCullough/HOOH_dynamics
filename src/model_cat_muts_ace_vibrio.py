'''

name:   model_cat_muts_ace_vibrio.py 

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
df_aceMHM = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'Cat_muts_MHMace', header = 0)
df_MHMnoC = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'Cat_muts_noC', header = 0)

df_all = df_aceMHM


#df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(hrs)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_a = df_all.loc[df_all['id'].str.contains('abiotic', case=False)].copy()  

df = df_all
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['abundance'] = df[['rep1','rep2','rep3']].mean(axis=1)
df['sigma'] = df[['rep1','rep2','rep3']].std(axis=1)

df['log_abundance'] = df[['log1','log2', 'log3']].mean(axis=1)
df['log_sigma'] = df[['log1','log2', 'log3']].std(axis=1)

df['log_sigma'] = 0.2
df.loc[df['organism'] == 'H', 'log_sigma'] = 0.08

##############
#setting working directory


#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[(df['treatment(nM)']==0)]  #assay 0 H 
df15 = df.loc[(df['treatment(nM)']==1500)]

df = df0

df0

## Reading in inits files for 0 and 400 models respectively
inits0 = pd.read_csv("../data/inits/MHMace_inits0.csv")


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




#setting how many MCMC chains you will run 
nits = 10000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params


# state variable names
snames = ['S','N'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.0002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.5})
#setting state variiable  prior guess
D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+5})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
D0_mean = inits0['D0'][0]
N0_mean = inits0['N0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_0H,
                          parameter_names=['k1','k2','D0','N0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          D0 = D0_prior.copy(),
                          N0  = N0_prior.copy(),
                          t_steps=1000,
                          D = D0_mean,
                          N = N0_mean,
                            )
    return M

def mono_0H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2 = params[0], params[1]
    D,N = max(y[0],0),max(y[1],0),
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dDdt = (k2 * N /( (ksp) + N) )*D     
    dNdt =  - (k2 * N /( (ksp) + N) )*D
    return [dDdt,dNdt]


def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)
#####################################
# Create and Run model on 0 and 400 df
#####################################

a0 = get_model(df0) 

posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod0 = a0.integrate()

a0res = get_residuals(a0)  #is this using the best fit or just a first run???


#########################################################
# graphing df and models together
#########################################################
c0 = 'deepskyblue'


# Set up graph for Dynamics and param histograms

fig1,ax1 = plt.subplots(1,3,figsize=[9,4]) #plot creation and config 
#set titles of subplots
fig1.suptitle('Vibrio in 0nM HOOH ') #full title config
fig1.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30) #shift white space for better fig view
ax1[0].set_xlabel('Time (days)')
ax1[1].set_title('Sh')
ax1[1].set_ylabel('Frequency')
ax1[1].set_xlabel('Parameter Value')
ax1[2].set_title('deltah')
ax1[2].set_ylabel('Frequency')
ax1[2].set_xlabel('Parameter Value (Logged)')


#plot dynamics of data and model for 0 assay 
ax1[0].plot(df.time,df.abundance, marker='o',color = c0, label = 'H data ') #data of 0 H assay
ax1[0].plot(mod0.time,mod0['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a0.plot_uncertainty(ax1[0],posteriors0,'H',100)

# plot histograms of params next to dynamics graphs
ax1[1].hist(((posteriors0.Sh)), facecolor=c0) #graphing Sh of 0 H assay 
ax1[2].hist(((posteriors0.deltah)), facecolor=c0) #graphing deltah of 0 H assay 

#config legends
l1 = ax1[0].legend(loc = 'lower right')
l1.draw_frame(False)


fig1.savefig('../figures/MHMace_vibrio_0_dynamics')



#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig2,ax2 = plt.subplots(1,2,sharex=True,figsize=[8,5]) #make plot
fig2.suptitle('Trace plots for Logged Params in 0H') #set main title 
fig2.subplots_adjust(right=0.90, wspace = 0.55, top = 0.90) #shift white space for better fig view
fig2.supxlabel('Model Iteration') #set overall x title 

ax2[0].set_ylabel('Log Sh')
ax2[1].set_ylabel('Log deltah')

#graphing iteration number vs parameter numbert logged 
ax2[0].scatter(posteriors0.iteration,np.log(posteriors0.Sh),color = c0)
ax2[1].scatter(posteriors0.iteration,np.log(posteriors0.deltah),color = c0)

#print out plot
fig2.savefig('../figures/MHMace_vibrio_0_TRACE')



#making and confing of residuals plot
fig3,ax3 = plt.subplots(1,2,figsize=[10,6])
fig3.suptitle('Model Dynamics and Output')



ax3[2].set_title('Residuals vs Model Fit Value  - abiotic 0')
ax3[2].set_ylabel('Model Value (H)')
ax3[2].set_xlabel('Residual')
#config legends for data differentialtion 
fig3.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)


#format fig  
ax3[0].set_title(S +' dynamics') #graph title for graph 1
ax3[0].semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('HOOH dynamics ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax3[1].set_xlabel('Time (days)') #settign x axis label for graph 1
ax3[0].set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax3[0].set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('HOOH (nM)')
#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax3[0].errorbar(df15[df15['organism']=='D']['time'],df15[df15['organism']=='D']['avg1'],yerr=df15[df15['organism']=='D']['std1'], marker='o', label = 'avg1')
ax3[0].errorbar(df15[df15['organism']=='D']['time'],df15[df15['organism']=='D']['avg2'],yerr=df15[df15['organism']=='D']['std2'], marker='v', label = 'avg2 ')
ax3[0].errorbar(df15[df15['organism']=='D']['time'],df15[df15['organism']=='D']['abundance'],yerr=df15[df15['organism']=='D']['sigma'], marker='d', label = 'MEAN ')

ax3[1].scatter(a0res['res'], a0res['abundance'],label = '0 H', color = c0) #where )



l3 = ax3.legend()
l3.draw_frame(False)
plt.legend()


#how to get residuals from all posterior runs not just best???

#print out plot
fig3.savefig('../figures/MHMace_vibrio_0_residuals')

plt.show()


pframe0 = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe0.to_csv("../data/inits/MHMace_inits0.csv")








##################################
#aVibrio in  spike df1500
##################################








inits15 = pd.read_csv("../data/inits/MHMace_inits1500.csv")


# state variable names
snames = ['D','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.00002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.6})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':12})
#setting state variiable  prior guess
D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':100})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
D0_mean = inits15['D0'][0]
N0_mean = inits15['N0'][0]
H0_mean = inits15['H0'][0]

# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 10000


#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','D0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
                          D0 = D0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          D = D0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    D,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dDdt = (k2 * N /( (ksp) + N) )*D - kdam*D*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*D
    dHdt = Sh- phi*D*H
    return [dDdt,dNdt,dHdt]





#####################################
# Create and Run model on 0 and 400 df
#####################################

a15 = get_model(df15) 

posteriors15 = a15.MCMC(chain_inits=inits15,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod15 = a15.integrate()

a15res = get_residuals(a15)  #is this using the best fit or just a first run???


#########################################################
# graphing df and models together
#########################################################
c15 = 'darkviolet'


# Set up graph for Dynamics and param histograms

fig4,ax4 = plt.subplots(1,3,figsize=[9,4]) #plot creation and config 
#set titles of subplots
fig4.suptitle('Cat Muts 1500 HOOH Model Output') #full title config
fig4.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30) #shift white space for better fig view
ax4[0].set_title('HOOH Dynamics')
ax4[0].set_ylabel('HOOH Concentration nM/mL')
ax4[0].set_xlabel('Time (days)')
ax4[1].set_title('Sh')
ax4[1].set_ylabel('Frequency')
ax4[1].set_xlabel('Parameter Value')
ax4[2].set_title('deltah')
ax4[2].set_ylabel('Frequency')
ax4[2].set_xlabel('Parameter Value (Logged)')


'''
#format fig  
ax3[0].set_title(S +' dynamics') #graph title for graph 1
ax3[0].semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('HOOH dynamics ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax3[1].set_xlabel('Time (days)') #settign x axis label for graph 1
ax3[0].set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax3[0].set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('HOOH (nM)')
#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df15[df15['organism']=='S']['time'],df15[df15['organism']=='S']['avg1'],yerr=df15[df15['organism']=='S']['std1'], marker='o', label = 'avg1')
ax0.errorbar(df15[df15['organism']=='S']['time'],df15[df15['organism']=='S']['avg2'],yerr=df15[df15['organism']=='S']['std2'], marker='v', label = 'avg2 ')
ax0.errorbar(df15[df15['organism']=='S']['time'],df15[df15['organism']=='S']['abundance'],yerr=df15[df15['organism']=='S']['sigma'], marker='d', label = 'MEAN ')
# graph 400 H assay even and odd avgs
ax1.errorbar(df15[df15['organism']=='H']['time'],df15[df15['organism']=='H']['avg1'],yerr=df15[df15['organism']=='H']['std1'], marker='o', label = 'avg1')
ax1.errorbar(df15[df15['organism']=='H']['time'],df15[df15['organism']=='H']['avg2'],yerr=df15[df15['organism']=='H']['std2'], marker='v', label = 'avg2')
ax1.errorbar(df15[df15['organism']=='H']['time'],df15[df15['organism']=='H']['abundance'],yerr=df15[df15['organism']=='H']['sigma'], marker='d', label = 'MEAN')

ax3[2].scatter(a0res['res'], a0res['abundance'],label = '0 H', color = c0) #where )


'''




#plot dynamics of data and model for 0 assay 
ax4[0].plot(df15.time,df15.abundance, marker='o',color = c15, label = 'H data ') #data of 0 H assay
ax4[0].plot(mod15.time,mod15['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a0.plot_uncertainty(ax4[0],posteriors15,'H',100)

# plot histograms of params next to dynamics graphs
ax4[1].hist(((posteriors15.Sh)), facecolor=c15) #graphing Sh of 0 H assay 
ax4[2].hist(((posteriors15.deltah)), facecolor=c15) #graphing deltah of 0 H assay 

#config legends
l4 = ax4[0].legend(loc = 'lower right')
l4.draw_frame(False)


fig4.savefig('../figures/MHMace_vibrio_1500_dynamics')


#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig5,ax5 = plt.subplots(1,2,sharex=True,figsize=[8,5]) #make plot
fig5.suptitle('Trace plots for Logged Params in 1500 nM H') #set main title 
fig5.subplots_adjust(right=0.90, wspace = 0.55, top = 0.90) #shift white space for better fig view
fig5.supxlabel('Model Iteration') #set overall x title 

ax5[0].set_ylabel('Log Sh')
ax5[1].set_ylabel('Log deltah')

#graphing iteration number vs parameter numbert logged 
ax5[0].scatter(posteriors15.iteration,np.log(posteriors15.Sh),color = c15)
ax5[1].scatter(posteriors15.iteration,np.log(posteriors15.deltah),color = c15)

#print out plot
fig5.savefig('../figures/MHMace_vibrio_1500_TRACE')


###############

#making and confing of residuals plot
fig6,ax6 = plt.subplots(figsize=[8,5])
fig6.suptitle('Residuals vs Model Fit Value  - abiotic 1500nM')
fig6.supylabel('Model Value (H)')
fig6.supxlabel('Residual')
#config legends for data differentialtion 
l6 = ax6.legend()
l6.draw_frame(False)

#plotting residual function output residual and abundance columns 
ax6.scatter(a15res['res'], a15res['abundance'],label = '1500nM H', color = c0) #where )

#how to get residuals from all posterior runs not just best???

#print out plot
fig6.savefig('../figures/MHMace_vibrio_1500_residuals')


###################


plt.show()


pframe15 = pd.DataFrame(a15.get_parameters(),columns=a15.get_pnames())
pframe15.to_csv("../data/inits/MHMace_inits1500.csv")

##############################
#vibrio modeling
##############################




# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print('\n Done my guy \n')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



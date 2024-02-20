
'''

name:   model_spiked_syn_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/HOOH_dynamics/src'
    
author: DKM

goal: Loop model and graph 0 H Pro assay and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop for 0 and 400 using different init files? (need H connected to 0 H first) 

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
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

df = df_mono
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

#strain slice trhough vol number selection 
vol = 52

#vol52 colors 
c0 = 'cornflowerblue'
c1 = 'darkorange'

#vol28 colors
#c0 = 'dodgerblue'
#c1 = 'tomato'

#vol53 colors WSyn CC9605 
#c0 = 'steelblue'
#c1 = 'chocolate'

#vol54 colors Syn WH7802 
#c0 = 'darkcyan'
#c1 = 'lightcoral'

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False) & (df['Vol_number']== vol)]  #assay 0 H 
df4 = df.loc[(df['assay'].str.contains('4', case=False)) & (df['Vol_number']== vol)]

df = df4

#####################################################
#config data in df from raw for odelib usefulness
#####################################################

#making avg columns of technical reps (std here only for graphing, not logged here)
#splitting df of Pro into 0 and 400 H assays 

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig2.suptitle('Syn Vol '+str(vol)+'  Monoculture in 400 nM HOOH')
fig2.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)

#format fig  
ax0.set_title('Syn Vol '+str(vol)+' dynamics') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('HOOH dynamics ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('HOOH (nM)')
#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df4[df4['organism']=='S']['time'],df4[df4['organism']=='S']['avg1'],yerr=df4[df4['organism']=='S']['std1'],color = 'b', marker='o', label = 'avg1')
ax0.errorbar(df4[df4['organism']=='S']['time'],df4[df4['organism']=='S']['avg2'],yerr=df4[df4['organism']=='S']['std2'], color = 'g',marker='v', label = 'avg2 ')
ax0.errorbar(df4[df4['organism']=='S']['time'],df4[df4['organism']=='S']['abundance'],yerr=df4[df4['organism']=='S']['sigma'], color = c0,marker='d', label = 'MEAN ')
# graph 400 H assay even and odd avgs
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['avg1'],yerr=df4[df4['organism']=='H']['std1'], color = 'b',marker='o', label = 'avg1')
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['avg2'],yerr=df4[df4['organism']=='H']['std2'], color = 'g',marker='v', label = 'avg2')
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['sigma'],color = c0, marker='d', label = 'MEAN')

plt.legend()

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits4 = pd.read_csv('../data/inits/syn_vol'+str(vol)+ '_inits4.csv')

#setting how many MCMC chains you will run 
nits = 10000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['S','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search

#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.00002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.6})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':12})
#setting state variiable  prior guess
S0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':100})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
S0_mean = inits4['S0'][0]
N0_mean = inits4['N0'][0]
H0_mean = inits4['H0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','S0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
                          S0 = S0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          S = S0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    S,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dSdt = (k2 * N /( (ksp) + N) )*S - kdam*S*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*S
    dHdt = 12- phi*S*H
    return [dSdt,dNdt,dHdt]

def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)

#df0.loc[:,'log_abundance'] = np.log(10**df0.log_abundance)

# get_models
a4 = get_model(df4) 

#broken here!!!!!!!!!!
# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,static_parameters =set(['k1','k2','N0']),print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()
a4res = get_residuals(a4)  #is this using the best fit or just a first run???

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################
###### fig set up

#########################################################
#graphing P model vs data and params histograms 
#########################################################

figall,axall = plt.subplots(2,3,figsize=[12,8])
figall.subplots_adjust(wspace=0.3,hspace=0.3)
ax4,ax5 = axall[0,:],axall[1,:]

# set up graph
#set titles and config graph 
ax4[1].set_title(r'$S_0$', fontsize = 14)
ax4[2].set_title(r'$\kappa_{dam}$', fontsize = 14)
ax4[0].semilogy()

ax4[0].set_ylabel('Cells (ml$^{-1}$)', fontsize = 14)
ax4[0].set_xlabel('Time (days)', fontsize = 14)
ax4[1].set_xlabel('Parameter Value', fontsize = 14)
ax4[1].set_ylabel('Frequency', fontsize = 14)
ax4[2].set_xlabel('Parameter Value', fontsize = 14)
ax4[2].set_ylabel('Frequency', fontsize = 14)

ax4[1].tick_params(axis='x', labelsize=14)
ax4[1].tick_params(axis='y', labelsize=14)
ax4[2].tick_params(axis='x', labelsize=14)
ax4[2].tick_params(axis='y', labelsize=14)

#shift fig subplots

#graph data, model, and uncertainty 
ax4[0].plot(df4[df4['organism']=='S']['time'], df4[df4['organism']=='S']['abundance'], color = c0, marker='o',label = r'$Synechococcus$')
ax4[0].plot(mod4.time,mod4['S'],color='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'S',100)

# plot histograms of parameter search results 
ax4[1].hist(posteriors4.S0, facecolor = c0)
ax4[2].hist(posteriors4.kdam, facecolor = c0)

#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
#show full graph 

#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 
ax5[1].set_title(r'$H_0$', fontsize = 14)
ax5[2].set_title('$\phi_{max}$', fontsize = 14)

ax5[0].semilogy()

ax5[0].set_ylabel(r'H$_2$O$_2$ concentration', fontsize = 14)
ax5[0].set_xlabel('Time (Days)', fontsize = 14)
ax5[1].set_xlabel('Parameter Value', fontsize = 14)
ax5[1].set_ylabel('Frequency', fontsize = 13)
ax5[2].set_xlabel('Parameter Value', fontsize = 14)
ax5[2].set_ylabel('Frequency', fontsize = 13)

ax5[1].tick_params(axis='x', labelsize=14)
ax5[1].tick_params(axis='y', labelsize=14)
ax5[2].tick_params(axis='x', labelsize=14)
ax5[2].tick_params(axis='y', labelsize=14)

#plot dynamics and models
ax5[0].plot(df4[df4['organism']=='H']['time'], df4[df4['organism']=='H']['abundance'], color = c1, marker='o',label = 'Hydrogen peroxide')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)

# plot histograms of parameter search results 
ax5[1].hist(posteriors4.H0, facecolor = c1)
ax5[2].hist(posteriors4.phi, facecolor = c1)

#make legends
l5 = ax5[0].legend(loc = 'upper right')
l5.draw_frame(False)
#show full graph 

for (ax,l) in zip(axall.flatten(),'abcdef'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

figall.savefig('../figures/syn_odelib4_Hparams')

#####residuals
fig6, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6)) #fig creationg of 1 by 2
fig6.suptitle('Syn in 400 H Model') #setting main title of fig

####### fig config and naming 

fig6.subplots_adjust(right=0.90, wspace = 0.45, left = 0.10, hspace = 0.20, bottom = 0.2)

ax0.semilogy()
ax0.set_title('Syn  dynamics ',fontsize = '16')
ax1.set_title('Model residuals',fontsize = '14')

ax0.set_ylabel('Data S value',fontsize = '14')
ax0.set_xlabel('Residual',fontsize = '14')

ax1.set_ylabel('Data H value',fontsize = '14')
ax1.set_xlabel('Residual',fontsize = '14')


ax0.scatter(a4res['res'], a4res['abundance'],label = '0H case')

ax1.scatter(a4res['res'], a4res['abundance'],label = '0H case')
#printing off graph
plt.show()



pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/syn_vol'+str(vol)+ '_inits4.csv')



# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')


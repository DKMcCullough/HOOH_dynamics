
'''

name:   model_spiked_pro_batch1.py 

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


######################################################
#reading in data and configureing 
#####################################################
df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'BCC_1-31-dataset', header = 1)

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

#####################################################
#config data in df from raw for odelib usefulness
#####################################################

#splitting df of Pro into 0 and 400 H assays 
 
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

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False) & (df['Vol_number']== 1)]  #assay 0 H 
df4 = df.loc[(df['assay'].str.contains('4', case=False)) & (df['Vol_number']== 1)]


df = df4

c0 = 'g'
c1 = 'darkgoldenrod'

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (9,5))
fig2.suptitle('Pro  Monoculture in 400 nM HOOH')
fig2.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)


#format fig  
ax0.set_title('Pro dynamics') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('HOOH dynamics ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('HOOH (nM)')

#graph 0 H assay even and odd avgs 
ax0.errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['avg1'],yerr=df4[df4['organism']=='P']['std1'], marker='o', label = 'avg1')
ax0.errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['avg2'],yerr=df4[df4['organism']=='P']['std2'], marker='v', label = 'avg2 ')
ax0.errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['abundance'],yerr=df4[df4['organism']=='P']['sigma'], marker='d', label = 'MEAN ')
# graph 400 H assay even and odd avgs
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['avg1'],yerr=df4[df4['organism']=='H']['std1'], marker='o', label = 'avg1')
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['avg2'],yerr=df4[df4['organism']=='H']['std2'], marker='v', label = 'avg2')
ax1.errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['sigma'], marker='d', label = 'MEAN')

plt.legend()

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits4 = pd.read_csv("../data/inits/pro_MIT9215_inits4_1.csv")

#setting how many MCMC chains you will run 
nits = 10000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['P','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.000002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.6})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':100})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
P0_mean = inits4['P0'][0]
N0_mean = inits4['N0'][0]
H0_mean = inits4['H0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','P0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
                          P0 = P0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          P = P0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    P,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P - kdam*P*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    dHdt = Sh- phi*P*H
    return [dPdt,dNdt,dHdt]


def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)

#####################################################
#modeling and fitting 
#####################################################

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
fig3, ax3 = plt.subplots(1,2,figsize = (9,5)) #fig creationg of 1 by 2
fig3.suptitle('Pro in 400 H Model') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.85, wspace = 0.50, hspace = 0.30)

ax3[0].semilogy()
ax3[1].semilogy()
ax3[0].set_title('Pro dynamics ', fontsize = 14)
ax3[1].set_title('HOOH dynamics', fontsize = 14)

ax3[0].set_ylabel('Cells (ml$^{-1}$)')
ax3[0].set_xlabel('Time (days)')
ax3[1].set_ylabel('HOOH (nM)')
ax3[1].set_xlabel('Time (days)')


#graphing data from df to see 2 different biological reps represented

ax3[0].errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['abundance'],yerr=df4[df4['organism']=='P']['std1'],c = c0, marker='o', label = 'Mean P')
ax3[0].plot(mod4.time,mod4['P'],color ='r',lw=1.5,label=' P model best fit')
a4.plot_uncertainty(ax3[0],posteriors4,'P',100)

ax3[1].errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['std1'],c = c1, marker='o', label = 'Mean H')
ax3[1].plot(mod4.time,mod4['H'],color ='r',lw=2.0,label=' H model best fit')
a4.plot_uncertainty(ax3[1],posteriors4,'H',100)

l3 = ax3[0].legend(loc = 'lower right')
l3.draw_frame(False)


#save graph

fig3.savefig('../figures/pro1_data_400')


#########################################################
#graphing P model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,3,figsize=[10,4])
#set titles and config graph 
fig4.suptitle('Pro Monoculture parameters in 400 HOOH ')
ax4[0].set_title('Pro dyanmics')
ax4[1].set_title('P0')
ax4[2].set_title('kdam')

ax4[0].semilogy()
ax4[0].set_ylabel('Cells (ml$^{-1}$)')
ax4[0].set_xlabel('Time (Days)')


ax4[1].set_xlabel('Parameter Value', fontsize = 12)
ax4[1].set_ylabel('Frequency', fontsize = 12)
ax4[2].set_xlabel('Parameter Value', fontsize = 12)
ax4[2].set_ylabel('Frequency', fontsize = 12)

#shift fig subplots
fig4.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)

ax4[0].set_ylim([100, 5000000])
#graph data, model, and uncertainty 
ax4[0].plot(df4[df4['organism']=='P']['time'], df4[df4['organism']=='P']['abundance'], color =  c0,marker='o',label = 'MEAN Pro')
ax4[0].plot(mod4.time,mod4['P'],color='r',lw=1.5,label=' Model P best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'P',100)

# plot histograms of parameter search results 
ax4[1].hist(posteriors4.P0, color =  c0)
ax4[2].hist(posteriors4.kdam,color = c0)

#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
#show full graph 
plt.show()
fig4.savefig('../figures/pro1_odelib4_Pparams')


#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 
fig5,ax5 = plt.subplots(1,3,figsize=[10,4])
fig5.suptitle('HOOH parmaters in 400 HOOH')
ax5[0].set_title('HOOH dynamics')
ax5[1].set_title('H0')
ax5[2].set_title('\u03C6')

ax5[0].set_ylabel('HOOH (nM)')
ax5[0].set_xlabel('Time (Days)')
fig5.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)

ax5[0].set_ylim([200, 500])
ax5[1].set_xlabel('Parameter Value', fontsize = 12)
ax5[1].set_ylabel('Frequency', fontsize = 12)
ax5[2].set_xlabel('Parameter Value', fontsize = 12)
ax5[2].set_ylabel('Frequency', fontsize = 12)

l5 = ax5[0].legend(loc = 'upper left')
l5.draw_frame(False)

ax5[0].plot(df4[df4['organism']=='H']['time'], df4[df4['organism']=='H']['abundance'], color =  c1,marker='o',label = 'H data ')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label=' Model H best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)


# plot histograms of parameter search results 
ax5[1].hist(posteriors4.H0,color =  c1)
ax5[2].hist(posteriors4.phi, color = c1)


#show full graph 

fig5.savefig('../figures/pro1_odelib4_Hparams')

##########################
#TRACE plot for death params
fig6,ax6 = plt.subplots(1,2,sharex=True,figsize=[8,4]) #make plot
fig6.suptitle('Trace plots for Params ', fontsize = 14) #set main title 
fig6.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view

ax6[0].set_title('kdam', fontsize = 14)
ax6[1].set_title('\u03C6', fontsize = 14)
ax6[0].set_ylabel('kdam value', fontsize = 12)
ax6[0].set_xlabel('Model iteration', fontsize = 12)
ax6[1].set_ylabel('\u03C6 value', fontsize = 12)
ax6[1].set_xlabel('Model iteration', fontsize = 12)
#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax6[0].scatter(posteriors4.iteration,posteriors4.kdam,color = c0)
ax6[1].scatter(posteriors4.iteration,posteriors4.phi,color = c1)


plt.show()
fig5.savefig('../figures/pro1_odelib4_Hparams')

#############
#residual graph 
#################

fig6, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6)) #fig creationg of 1 by 2
fig6.suptitle('Pro in 400 H Model') #setting main title of fig

####### fig config and naming 

fig6.subplots_adjust(right=0.90, wspace = 0.45, left = 0.10, hspace = 0.20, bottom = 0.2)

ax0.semilogy()
ax0.set_title('Pro  dynamics ',fontsize = '16')
ax1.set_title('Model residuals',fontsize = '14')

ax0.set_ylabel('Data P value',fontsize = '14')
ax0.set_xlabel('Residual',fontsize = '14')

ax1.set_ylabel('Data H value',fontsize = '14')
ax1.set_xlabel('Residual',fontsize = '14')


ax0.scatter(a4res['res'], a4res['abundance'],label = '0H case')

ax1.scatter(a4res['res'], a4res['abundance'],label = '0H case')
#printing off graph
plt.show()


#save over best params to new inits
pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/pro_MIT9215_inits4_1.csv')


# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')

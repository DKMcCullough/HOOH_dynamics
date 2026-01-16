
'''

name:   model_spiked_EZ55.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/HOOH_dynamics/src'
    
author: DKM

goal: Loop model and graph 0 H Pro assay and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop for 0 and 400 using different init files? (need H connected to 0 H first) 

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
df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'EZ55_REDO', header = 0)
df_all.rename(columns = {'Time (hours)':'time'},inplace = True) 

dfw = df_all
dfw['time'] = dfw['time']      #changing time foorm hours to days
dfw = dfw[dfw.time >0.3] #to keep growth bump in last days to thro off death (kdma) range.

dfw['log1'] = np.log(dfw['rep1'])
dfw['log2'] = np.log(dfw['rep2'])
dfw['log3'] = np.log(dfw['rep3'])
#dfw['log4'] = np.log(dfw['rep4'])
#dfw['log5'] = np.log(dfw['rep5'])
#dfw['log6'] = np.log(dfw['rep6'])
dfw['abundance'] = dfw[['rep1','rep2','rep3']].mean(axis=1)
dfw['sigma'] = dfw[['rep1','rep2','rep3']].std(axis=1)
dfw['log_abundance'] = dfw[['log1','log2', 'log3']].mean(axis=1)
dfw['log_sigma'] = dfw[['log1','log2', 'log3']].std(axis=1)

dfw['log_sigma'] = 0.1
dfw.loc[dfw['organism'] == 'H', 'log_sigma'] = 0.08
dfw.time = dfw.time - 0.3

dfa =dfw.loc[(dfw['strain'].str.contains('abio', case=False))]
dfb =dfw.loc[(~dfw['strain'].str.contains('abio', case=False))]

df_all.rename(columns = {'Time (days)':'time'},inplace = True) 

df_all.rename(columns = {'Time (days)':'time'},inplace = True) 


c0 = 'b'
c1 = 'orange'

#Het_55_inits4.csv

IDs = dfb.ID.unique()

for i in IDs:
    df =dfb[dfb['ID']== i]

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
    inits4 = pd.read_csv('../data/inits/Het_55_inits_'+str(i)+'.csv')
    inits0 = pd.read_csv('../data/inits/Het_55_inits_'+str(IDs[0])+'.csv')

#setting how many MCMC chains you will run 
    nits = 100000
    # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
    snames = ['D','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
    pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
    k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.000002})
    k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':inits0['k2'][0]})
    kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.00002})
    phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.006})
    Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
    D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
    N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
    H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':300})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
    D0_mean = inits4['D0'][0]
    N0_mean = inits4['N0'][0]
    H0_mean = inits4['H0'][0]

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
        dHdt = 8 - phi*D*H
        return [dDdt,dNdt,dHdt]

#####################################################
#modeling and fitting 
#####################################################

# get_models
    a4 = get_model(df) 

#broken here!!!!!!!!!!
# do fitting
    posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,static_parameters =set(['k1','k2','N0']),print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
    mod4 = a4.integrate()

#########################################################
#graphing P model vs data and params histograms 
#########################################################

# set up graph
    figall,axall = plt.subplots(2,3,figsize=[12,8])
    figall.subplots_adjust(wspace=0.3,hspace=0.3)
    ax4,ax5 = axall[0,:],axall[1,:]
#set titles and config graph 
    ax4[0].semilogy()

    ax4[0].set_ylabel('Cells (mL$^{-1}$)', fontsize = 12)
    ax4[0].set_xlabel('Time (days)', fontsize = 12)
    ax4[1].set_xlabel('log$_{10}$ Initial cell density \n ($P_{i,0}$, cells mL$^{-1}$)', fontsize = 12)
    ax4[1].set_ylabel('Probability density', fontsize = 12)
    ax4[2].set_xlabel('log$_{10}$ Damage rate \n ($\kappa_{dam,i}$, mL pmol$^{-1}$ day$^{-1}$)', fontsize = 12)
    ax4[2].set_ylabel('Probability density', fontsize = 12)

    #ax4[0].set_ylim([100, 5000000])
#graph data, model, and uncertainty 
    ax4[0].plot(df[df['organism']=='D']['time'], df[df['organism']=='D']['abundance'], color =  c0,marker='o',label = r'$Alteromonas$ $EZ55$')
    ax4[0].plot(mod4.time,mod4['D'],color='r',lw=1.5,label='Model best fit')
    a4.plot_uncertainty(ax4[0],posteriors4,'D',100)

# plot histograms of parameter search results 
    hp.sns.kdeplot(np.log10(posteriors4.D0),color=c0,fill=True,ax=ax4[1],linewidth=3)
    hp.sns.kdeplot(np.log10(posteriors4.kdam),color=c0,fill=True,ax=ax4[2],linewidth=3)

#make legends
    l4 = ax4[0].legend(loc = 'lower left')
    l4.draw_frame(False)
#show full graph 

    for (ax,l) in zip(axall.flatten(),'abcdef'):
        ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 

    ax5[0].semilogy()

    ax5[0].set_ylabel(r'H$_2$O$_2$ concentration (pmol mL$^{-1}$)', fontsize=12)
    ax5[0].set_xlabel('Time (Days)', fontsize = 12)
    ax5[1].set_xlabel('log$_{10}$ Initial H$_2$O$_2$ concentration \n ($H_0$, pmol mL$^{-1}$)', fontsize = 12)
    ax5[1].set_ylabel('Probability density', fontsize = 12)
    ax5[2].set_xlabel('log$_{10}$ Detoxification rate \n ($\phi_{det,i}$, pmol cell$^{-1}$ day$^{-1}$)', fontsize = 12)
    ax5[2].set_ylabel('Probability density', fontsize = 12)

    ax5[0].plot(df[df['organism']=='H']['time'], df[df['organism']=='H']['abundance'], color =  c1,marker='o',label = 'Hydrogen peroxide')
    ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label='Model best fit')
    a4.plot_uncertainty(ax5[0],posteriors4,'H',100)

    l5 = ax5[0].legend(loc = 'lower left')
    l5.draw_frame(False)

# plot histograms of parameter search results 
    hp.sns.kdeplot(np.log10(posteriors4.H0),color=c1,fill=True,ax=ax5[1],linewidth=3)
    hp.sns.kdeplot(np.log10(posteriors4.phi),color=c1,fill=True,ax=ax5[2],linewidth=3)

#show full graph 
    figall.savefig('../figures/EZ55_h2o2_'+str(i))

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


#save over best params to new inits
    pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
    pframe.to_csv('../data/inits/Het_55_inits_'+str(i)+'.csv')

    posteriors4.to_csv('../data/posteriors/pos_EZ55_'+str(i)+ '.csv')

# 'program finished' flag
#plt.show()


'''

name:   model_spiked_pro_cocultre.py 

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
df_P = df_all[df_all['organism']== 'P']

df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  

#####################################################
#config data in df from raw for odelib usefulness
#####################################################

dfw = df_co

dfw['log1'] = np.log(dfw['rep1'])
dfw['log2'] = np.log(dfw['rep2'])
dfw['log3'] = np.log(dfw['rep3'])
dfw['log4'] = np.log(dfw['rep4'])
dfw['avg1'] = dfw[['rep1', 'rep3']].mean(axis=1)
dfw['avg2'] = dfw[['rep2', 'rep4']].mean(axis=1)
dfw['abundance'] = dfw[['rep1','rep2','rep3', 'rep4']].mean(axis=1)
dfw['std1'] = dfw[['rep1', 'rep3']].std(axis=1)
dfw['std2'] = dfw[['rep2', 'rep4']].std(axis=1)
dfw['sigma'] = dfw[['rep1','rep2','rep3', 'rep4']].std(axis=1)

dfw['lavg1'] = dfw[['log1', 'log3']].mean(axis=1) #making logged avg columns in df for odelib to have log_abundance to use for posterior calcs
dfw['lavg2'] = dfw[['log2', 'log4']].mean(axis=1)
dfw['log_abundance'] = dfw[['log1','log2', 'log3','log4']].mean(axis=1)
dfw['stdlog1'] = dfw[['log1', 'log3']].std(axis=1) #taking stdv of logged reps
dfw['stdlog2'] = dfw[['log2', 'log4']].std(axis=1)
dfw['log_sigma'] = dfw[['log1','log2', 'log3','log4']].std(axis=1)

dfw['log_sigma'] = 0.2
dfw.loc[dfw['organism'] == 'H', 'log_sigma'] = 0.08


###### Slicing by coculture (all different Syn Vol numbers)


#coculture #s coculture_1,28     1,52    1,53         1,54       1,55      

A = ('coculture_1,54')

df = dfw[dfw['assay']== A]

df[['log_abundance']] = df[['log_abundance']].astype(float)


###########
#colors for graphing
############3

c0 = 'c'
c1 = 'goldenrod'
c2 = 'dodgerblue'



#####################################################
#   model param and state variable set up 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits4 = pd.read_csv('../data/inits/'+str(A)+'_pro_inits4.csv')

#setting how many MCMC chains you will run 
nits = 1000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['P','S','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.000002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.6})
kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
kdams_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
phis_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.06})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':350})
S0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})

#still not sure what part of fitting algor this is used for
P0_mean = inits4['P0'][0]
S0_mean = inits4['S0'][0]
N0_mean = inits4['N0'][0]
H0_mean = inits4['H0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                          parameter_names=['k1','k2','kdam','phi','Sh','kdams','phis','P0','N0','H0','S0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          kdam = kdam_prior.copy(),
                          phi = phi_prior.copy(),
                          Sh = Sh_prior.copy(),
                          kdams = kdams_prior.copy(),
                          phis = phis_prior.copy(),
                          P0 = P0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          S0 = S0_prior.copy(),
                          t_steps=1000,
                          P = P0_mean,
                          S = S0_mean,
                          N = N0_mean,
                          H = H0_mean,
                            )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh,kdams,phis = params[0], params[1], params[2],params[3], params[4],params[5],params[6]
    P,S,N,H  = max(y[0],0),max(y[1],0),max(y[2],0),max(y[3],0)
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P - kdam*P*H    
    dSdt = (k2 * N /( (ksp) + N) )*S - kdams*S*H  
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    dHdt = Sh- phi*P*H - phis*S*H
    return [dPdt,dSdt, dNdt,dHdt]


def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)

#####################################################
#modeling and fitting 
#####################################################

# get_models
a4 = get_model(df) 


# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,static_parameters =set(['k1','k2','N0']),print_report=True) #, )

# run model with optimal params
mod4 = a4.integrate()


#####################################################
# graphing model vs data 
#####################################################

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1,ax2)= plt.subplots(1,3,figsize = (10,5))
fig2.suptitle(' '+str(A)+' in 400 nM HOOH')
fig2.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)


#format fig  
ax0.set_title('Pro dynamics') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('HOOH dynamics ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax2.set_title('Pro dynamics') #graph title for graph 1
ax2.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('HOOH (nM)')
ax2.set_xlabel('Time (days)') #settign x axis label for graph 1
ax2.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 

#graph P
ax0.errorbar(df[df['organism']=='P']['time'],df[df['organism']=='P']['avg1'],yerr=df[df['organism']=='P']['std1'], marker='o', label = 'avg1')
ax0.errorbar(df[df['organism']=='P']['time'],df[df['organism']=='P']['avg2'],yerr=df[df['organism']=='P']['std2'], marker='v', label = 'avg2 ')
ax0.errorbar(df[df['organism']=='P']['time'],df[df['organism']=='P']['abundance'],yerr=df[df['organism']=='P']['sigma'], marker='d', label = 'MEAN P')
# graph H
ax1.errorbar(df[df['organism']=='H']['time'],df[df['organism']=='H']['avg1'],yerr=df[df['organism']=='H']['std1'], marker='o', label = 'avg1')
ax1.errorbar(df[df['organism']=='H']['time'],df[df['organism']=='H']['avg2'],yerr=df[df['organism']=='H']['std2'], marker='v', label = 'avg2')
ax1.errorbar(df[df['organism']=='H']['time'],df[df['organism']=='H']['abundance'],yerr=df[df['organism']=='H']['sigma'], marker='d', label = 'MEAN H')
#graph S 
ax2.errorbar(df[df['organism']=='S']['time'],df[df['organism']=='S']['avg1'],yerr=df[df['organism']=='S']['std1'], marker='o', label = 'avg1')
ax2.errorbar(df[df['organism']=='S']['time'],df[df['organism']=='S']['avg2'],yerr=df[df['organism']=='S']['std2'], marker='v', label = 'avg2')
ax2.errorbar(df[df['organism']=='S']['time'],df[df['organism']=='S']['abundance'],yerr=df[df['organism']=='S']['sigma'], marker='d', label = 'MEAN S')

l2 = ax0.legend(loc = 'upper right')
l2.draw_frame(False)


plt.legend()


########################
#Show model vs mean data 
########################

###### fig set up
fig3, ax3 = plt.subplots(1,3,figsize = (10,4)) #fig creationg of 1 by 2
fig3.suptitle('Pro '+str(A)+' in 400 H Model') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)

ax3[0].semilogy()
ax3[1].semilogy()
ax3[2].semilogy()

ax3[0].set_title('Pro dynamics ', fontsize = 14)
ax3[1].set_title('HOOH dynamics', fontsize = 14)
ax3[2].set_title('Syn dynamics', fontsize = 14)

ax3[0].set_ylabel('Cells (ml$^{-1}$)')
ax3[0].set_xlabel('Time (days)')
ax3[1].set_ylabel('HOOH (nM)')
ax3[1].set_xlabel('Time (days)')
ax3[2].set_ylabel('Cells (ml$^{-1}$)')
ax3[2].set_xlabel('Time (days)')

#graphing data from df to see 2 different biological reps represented

ax3[0].errorbar(df[df['organism']=='P']['time'],df[df['organism']=='P']['abundance'],yerr=df[df['organism']=='P']['std1'],c = c0, marker='o', label = 'Data Mean')
ax3[0].plot(mod4.time,mod4['P'],color ='r',lw=1.5,label=' Model best fit')
a4.plot_uncertainty(ax3[0],posteriors4,'P',100)

ax3[1].errorbar(df[df['organism']=='H']['time'],df[df['organism']=='H']['abundance'],yerr=df[df['organism']=='H']['std1'],c = c1, marker='o', label = 'Mean H')
ax3[1].plot(mod4.time,mod4['H'],color ='r',lw=2.0,label=' H model best fit')
a4.plot_uncertainty(ax3[1],posteriors4,'H',100)

ax3[2].errorbar(df[df['organism']=='S']['time'],df[df['organism']=='S']['abundance'],yerr=df[df['organism']=='S']['std1'],c = c2, marker='o', label = 'Mean S')
ax3[2].plot(mod4.time,mod4['S'],color ='r',lw=2.0,label=' S model best fit')
a4.plot_uncertainty(ax3[2],posteriors4,'S',100)

l3 = ax3[0].legend(loc = 'upper right')
l3.draw_frame(False)


#save graph

fig3.savefig('../figures/pro_coculture_'+str(A)+'data_')


#########################################################
#graphing P model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,4,figsize=[10,4])
#set titles and config graph 
fig4.suptitle(' '+str(A)+' Pro params ')
ax4[0].set_title('P0')
ax4[1].set_title('kdam')
ax4[2].set_title('\u03C6')
ax4[2].set_title('H0')

ax4[0].set_xlabel('P0 Value', fontsize = 12)
ax4[0].set_ylabel('Frequency', fontsize = 12)
ax4[1].set_xlabel('Kdam Value', fontsize = 12)
ax4[1].set_ylabel('Frequency', fontsize = 12)
ax4[2].set_xlabel('\u03C6 Value', fontsize = 12)
ax4[2].set_ylabel('Frequency', fontsize = 12)
ax4[3].set_xlabel('H0 Value', fontsize = 12)
ax4[3].set_ylabel('Frequency', fontsize = 12)
#shift fig subplots
fig4.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)


# plot histograms of parameter search results 
ax4[0].hist(posteriors4.P0, color =  c0)
ax4[1].hist(posteriors4.kdam,color = c0)
ax4[2].hist(posteriors4.phi,color = c1)
ax4[3].hist(posteriors4.H0,color = c1)

#show full graph 
plt.show()
fig4.savefig('../figures/pro_'+str(A)+'_odelib4_Pparams')


#########################################################
#graphing S params
#########################################################

# set up graph
fig5,ax5 = plt.subplots(1,4,figsize=[10,4])
#set titles and config graph 
fig5.suptitle(' '+str(A)+' Syn params ')
ax5[0].set_title('S0')
ax5[1].set_title('kdams')
ax5[2].set_title('\u03C6s')
ax5[3].set_title('H0')

ax5[0].set_xlabel('S0 Value', fontsize = 12)
ax5[0].set_ylabel('Frequency', fontsize = 12)
ax5[1].set_xlabel('Kdams Value', fontsize = 12)
ax5[1].set_ylabel('Frequency', fontsize = 12)
ax5[2].set_xlabel('\u03C6s Value', fontsize = 12)
ax5[2].set_ylabel('Frequency', fontsize = 12)
ax5[3].set_xlabel('H0 Value', fontsize = 12)
ax5[3].set_ylabel('Frequency', fontsize = 12)
#shift fig subplots
fig5.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)


# plot histograms of parameter search results 
ax5[0].hist(posteriors4.S0, color =  c2)
ax5[1].hist(posteriors4.kdams,color = c2)
ax5[2].hist(posteriors4.phis,color = c2)
ax5[3].hist(posteriors4.H0,color = c1)

#show full graph 
plt.show()
fig5.savefig('../figures/pro_'+str(A)+'_odelib4_Synparams')



##########################
#TRACE plot for death params
fig6,ax6 = plt.subplots(2,2,sharex=True,figsize=[9,7]) #make plot
fig6.suptitle('Trace plots for '+str(A)+' Params ', fontsize = 14) #set main title 
fig6.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, wspace=0.45, hspace=0.25) #shift white space for better fig view

ax6[0,0].set_title('kdam', fontsize = 14)
ax6[0,1].set_title('\u03C6', fontsize = 14)
ax6[1,0].set_title('kdams', fontsize = 14)
ax6[1,1].set_title('\u03C6s', fontsize = 14)
ax6[0,0].set_ylabel('kdam value', fontsize = 12)
ax6[0,0].set_xlabel('Model iteration', fontsize = 12)
ax6[0,1].set_ylabel('\u03C6 value', fontsize = 12)
ax6[0,1].set_xlabel('Model iteration', fontsize = 12)
ax6[1,0].set_ylabel('kdams value', fontsize = 12)
ax6[1,0].set_xlabel('Model iteration', fontsize = 12)
ax6[1,1].set_ylabel('\u03C6s value', fontsize = 12)
ax6[1,1].set_xlabel('Model iteration', fontsize = 12)
#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax6[0,0].scatter(posteriors4.iteration,posteriors4.kdam,color = c0)
ax6[0,1].scatter(posteriors4.iteration,posteriors4.phi,color = c1)
ax6[1,0].scatter(posteriors4.iteration,posteriors4.kdams,color = c0)
ax6[1,1].scatter(posteriors4.iteration,posteriors4.phis,color = c1)


plt.show()
fig6.savefig('../figures/'+str(A)+'_odelib4_trace')



#save over best params to new inits
pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/'+str(A)+'_pro_inits4.csv')


# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')

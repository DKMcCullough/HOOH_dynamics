'''

name:   model_vibrio_ace.py

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


#df0 = df[df['id'] =='Growth_MHMacetate_0nMspike'] #set id of assay we are fitting
df0 = df[df['id'] == 'Stationary_MHMacetate_0nMspike']
#df = df[df['id'] =='CatMutsNoC_stationary'] #set id of assay we are fitting

  #assay 0 H 
df15 = df[df['id'] =='Stationary_MHMacetate_1500nMspike'] #assay of spiked HOOH


#st  = 'WT'
st  = 'WT'


df = df0[df0['strain']==st]



inits0 = pd.read_csv("../data/inits/MHMace_inits0.csv")


    #####################################################
    #functions  for modeling and graphing model uncertainty 
    #####################################################

pw = 1

    #setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':.200})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.5})
    #setting state variiable  prior guess
D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+4})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+5})
    #pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

    #still not sure what part of fitting algor this is used for
D0_mean = inits0['D0'][0]
N0_mean = inits0['N0'][0]
snames = ['D','N']

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
    #print('D, dDdt, N, dNdt')
    #print(D, dDdt, N, dNdt)
    return [dDdt,dNdt]

    
    #####################################################
    
    #find closest time 
def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)
    
    
    
    # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 1000
    
    
    #####################################
    # Create and Run model on 0 and 400 df
    #####################################
    
a0 = get_model(df) 
    
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
    
    # run model with optimal params
mod0 = a0.integrate()
    
a0res = get_residuals(a0)  #is this using the best fit or just a first run???
    
    
    #########################################################
    # graphing df and models together
    #########################################################
c0 = 'c'
    
    
    # Set up graph for Dynamics and param histograms
    
fig1, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6)) #fig creationg of 1 by 2
fig1.suptitle(str(st)+' in 0 H Model') #setting main title of fig
    
    ####### fig config and naming 
    
fig1.subplots_adjust(right=0.9, wspace = 0.45, hspace = 0.20)
    
ax0.semilogy()
ax0.set_title(str(st) +' dynamics ')
ax1.set_title('Model residuals')
    
ax0.set_xlabel('Time (hrs)')
ax0.set_ylabel('Cells (ml-1)')
ax1.set_ylabel('Data D value')
ax1.set_xlabel('Residual')


    #graphing data from df to see 2 different biological reps represented
    
ax0.errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['rep1'],color = 'lightcoral', marker='^', label = 'rep1')
ax0.errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['rep2'],color = 'firebrick' ,marker='.', label = 'rep2')
ax0.plot(df[df['organism']== 'D']['time'], df[df['organism']== 'D']['rep3'], color = 'lightsalmon', marker='o',label = 'rep3')
ax0.errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['abundance'],yerr=df[df['organism']== 'D']['sigma'],color = c0, marker='d', label = 'MEAN '+str(st))

ax0.plot(mod0.time,mod0['D'],c='r',lw=1.5,label=' model best fit')

a0.plot_uncertainty(ax0,posteriors0,'D',100)

ax1.scatter(a0res['res'], a0res['abundance'],color = c0,label = '0H case')

l1 = ax0.legend(loc = 'lower right')
l1.draw_frame(False)
    
    
    
fig1.savefig('../figures/MHMace_WT_stat_dynamics')
    
    
    
    ###################
fig2,ax2 = plt.subplots(1,3, figsize=[8,5]) 
fig2.suptitle(('MHMace Model Dynamics for '+str(st) )) #set main title 
ax2[0].set_title((str(st) +' dynamics'), fontsize = 12)
ax2[1].set_title('D0', fontsize = 12)
ax2[2].set_title('k2', fontsize = 12)
ax2[0].semilogy()
#ax2[0].set_xlim(xmin =0 , xmax =35)
ax2[1].set_xlabel('Parameter Value', fontsize = 12)
ax2[1].set_ylabel('Frequency', fontsize = 12)
ax2[1].xaxis.set_label_coords(0.88, -0.2)
ax2[0].set_xlabel('Time (days)', fontsize = 14)
    #make legends

    #shift fig subplots
fig2.subplots_adjust(right=0.95, wspace = 0.45, left = 0.05, hspace = 0.20, bottom = 0.2)
    
    
    #graph data, model, and uncertainty 
#ax2[0].plot(df[df['organism']== 'D']['time'], df[df['organism']== 'D']['abundance'], marker='>',label = str(s) +' MEAN in 0 H ')
ax2[0].errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['abundance'],yerr=df[df['organism']== 'D']['sigma'], color = c0,marker='^', label = 'MEAN WT')
ax2[0].plot(mod0.time,mod0['D'],c='r',lw=1.5,label=' Model D best fit')
a0.plot_uncertainty(ax2[0],posteriors0,'D',100)
    
    
    # plot histograms of parameter search results 
ax2[1].hist(posteriors0.D0,color = c0)
ax2[1].tick_params(axis='x', labelsize=14)
ax2[1].tick_params(axis='y', labelsize=14)
    
ax2[2].hist(posteriors0.k2,color = c0)
ax2[2].tick_params(axis='x', labelsize=14)
ax2[2].tick_params(axis='y', labelsize=14)

l2 = ax2[0].legend(loc = 'upper left')
l2.draw_frame(False)

fig2.savefig('../figures/'+str(st) +'_MHMace_0H_odelib')

    

#################################
#graphing trace of Param values
##################################
#crating and config of fig 3
fig3,ax3 = plt.subplots(1,2, figsize=[8,5]) #make plot
fig3.suptitle((str(st) +' MHMnoC Trace plots in 0H')) #set main title 
fig3.subplots_adjust(right=0.90, wspace = 0.55, top = 0.90) #shift white space for better fig view
fig3.supxlabel('Model Iteration') #set overall x title 

ax3[0].set_ylabel('D0')
ax3[1].set_ylabel('k2')

#graphing iteration number vs parameter numbert logged 
ax3[0].scatter(posteriors0.iteration,(posteriors0.D0),color = c0)
ax3[1].scatter(posteriors0.iteration,(posteriors0.k2),color = c0)

#print out plot
fig3.savefig('../figures/MHMace_'+str(st) +'_0_TRACE')


plt.show()


pframe0 = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe0.to_csv("../data/inits/MHMace_inits0.csv")

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################


##################################
#Vibrio in  spiked df1500
##################################

df = df15[df15['strain']==st]

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


# get_models
a4 = get_model(df) 

#broken here!!!!!!!!!!
# do fitting
posteriors4 = a4.MCMC(chain_inits=inits15,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )# static_parameters =set(['k1','k2','N0'])
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()


#####################################################
# graphing model vs data in 0 H and associated error
#####################################################

###### fig set up
fig3, ax3 = plt.subplots(1,2,figsize = (8,5)) #fig creationg of 1 by 2
fig3.suptitle('Vibrio in 400 H Model MHMace') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.85, wspace = 0.40, hspace = 0.30)

ax3[0].semilogy()
ax3[1].semilogy()
ax3[0].set_title(str(st) +' vibrio dynamics ')
ax3[1].set_title('HOOH dynamics')

ax3[0].set_xlabel('days')
ax3[0].set_ylabel('cell concentration')
ax3[1].set_ylabel('HOOH concentration')
ax3[1].set_xlabel('days ')

l3 = ax3[0].legend(loc = 'lower right')
l3.draw_frame(False)

ax3[0].semilogy()
ax3[1].semilogy()
#graphing data from df to see 2 different biological reps represented

ax3[0].errorbar(df[df['organism']=='D']['time'],df[df['organism']=='D']['abundance'],yerr=df[df['organism']=='D']['sigma'], marker='o', label = 'Vibrio MEAN')
ax3[0].plot(mod4.time,mod4['D'],color ='r',lw=1.5,label=' D model best fit')
a4.plot_uncertainty(ax3[0],posteriors4,'D',100)

ax3[1].errorbar(df[df['organism']=='H']['time'],df[df['organism']=='H']['abundance'],yerr=df[df['organism']=='H']['sigma'], marker='o', label = 'Mean H')
ax3[1].plot(mod4.time,mod4['H'],color ='r',lw=2.0,label=' H model best fit')
a4.plot_uncertainty(ax3[1],posteriors4,'H',100)

#ax1.scatter(a0res['res'], a0res['abundance'],label = '0H case')
#printing off graph
plt.show()




#########################################################
#graphing P model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,3,figsize=[10,5])
#set titles and config graph 
fig4.suptitle('Vibrio Monoculture in 400 HOOH MHMace', fontsize = 16)
ax4[0].set_title(str(st) +' Dynamics', fontsize = 14)
ax4[1].set_title('D0', fontsize = 14)
ax4[2].set_title('kdam', fontsize = 14)

ax4[0].semilogy()

fig4.subplots_adjust(right=0.80, wspace = 0.25, hspace = 0.20)

ax4[1].set_xlabel('Parameter Value', fontsize = 14)
ax4[1].set_ylabel('Frequency', fontsize = 14)

ax4[0].set_ylabel('Cells (ml-1)')
ax4[0].set_xlabel('Time (days)', fontsize = 14)
#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)
ax4[1].tick_params(axis='x', labelsize=14)
ax4[1].tick_params(axis='y', labelsize=14)
ax4[2].tick_params(axis='x', labelsize=14)
ax4[2].tick_params(axis='y', labelsize=14)

#shift fig subplots


fig4.subplots_adjust(right=0.95, wspace = 0.45, left = 0.05, hspace = 0.30, bottom = 0.2)

#graph data, model, and uncertainty 
ax4[0].errorbar(df[df['organism']=='D']['time'],df[df['organism']=='D']['abundance'],yerr=df[df['organism']=='D']['sigma'], marker='o', label = 'Vibrio MEAN')
ax4[0].plot(mod4.time,mod4['D'],color='r',lw=1.5,label=' Model D best fit')
a4.plot_uncertainty(ax4[0],posteriors4,'D',100)

# plot histograms of parameter search results 
ax4[1].hist(posteriors4.D0)
ax4[2].hist(posteriors4.kdam)

#show full graph 
plt.show()
fig4.savefig('../figures/syn_odelib4_Sparams')


#########################################################
#graphing H model vs data and params histograms 
#########################################################

#HOOH dynamics 
fig5,ax5 = plt.subplots(1,3,figsize=[10,5])
fig5.suptitle('HOOH parmaters ', fontsize = 14)
ax5[0].set_title('HOOH Dynamics', fontsize = 14)
ax5[1].set_title('H0', fontsize = 14)
ax5[2].set_title('phi', fontsize = 14)

ax5[0].semilogy()
fig5.subplots_adjust(right=0.85, wspace = 0.45, hspace = 0.20)
ax5[1].set_xlabel('Parameter Value', fontsize = 14)
ax5[1].set_ylabel('Frequency', fontsize = 13)
ax5[0].set_ylabel('HOOH concentration', fontsize = 14)
ax5[0].set_xlabel('Time (Days)', fontsize = 14)
#make legends
l5 = ax5[0].legend(loc = 'upper left')
l5.draw_frame(False)
ax5[1].tick_params(axis='x', labelsize=14)
ax5[1].tick_params(axis='y', labelsize=14)
ax5[2].tick_params(axis='x', labelsize=14)
ax5[2].tick_params(axis='y', labelsize=14)

ax5[0].errorbar(df[df['organism']=='H']['time'],df[df['organism']=='H']['abundance'],yerr=df[df['organism']=='H']['sigma'], marker='o', label = 'HOOH MEAN')
ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label=' Model H best fit')
a4.plot_uncertainty(ax5[0],posteriors4,'H',100)


# plot histograms of parameter search results 
ax5[1].hist(posteriors4.H0)
ax5[2].hist(posteriors4.phi)


#show full graph 
plt.show()


pframe15 = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe15.to_csv("../data/inits/MHMace_inits1500.csv")


plt.show()

# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print('\n Done my guy \n')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')




'''

name:   model_syn_batch.py 

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

plt.rcParams["font.family"] = "Times New Roman"

######################################################
#reading in data and configureing 
#####################################################

df = hp.get_data('coculture')

vol = int(sys.argv[1])
if vol == 52:
    #vol52 colors Syn WH7803
    c0 = 'cornflowerblue'
    c1 = 'darkorange'
elif vol == 28:
    #vol28 colors Syn WH7803 also 
    c0 = 'dodgerblue'
    c1 = 'tomato'
elif vol == 53:
    #vol53 colors WSyn CC9605 
    c0 = 'steelblue'
    c1 = 'chocolate'
elif vol == 54:
    #vol54 colors Syn WH7802 
    c0 = 'darkcyan'
    c1 = 'lightcoral'

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False) & (df['Vol_number']== vol)]  #assay 0 H 
df4 = df.loc[(df['assay'].str.contains('4', case=False)) & (df['Vol_number']== vol)]

# quantify uncertainty
sigma0H = hp.get_uncertainty(df0[df0.organism=='H'])
sigma0S = hp.get_uncertainty(df0[df0.organism=='S'])
sigma4H = hp.get_uncertainty(df4[df4.organism=='H'])
sigma4S = hp.get_uncertainty(df4[df4.organism=='S'])
df0.loc[df['organism'] == 'H', 'log_sigma'] = sigma0H
df0.loc[df['organism'] == 'S', 'log_sigma'] = sigma0S
df4.loc[df['organism'] == 'H', 'log_sigma'] = sigma4H
df4.loc[df['organism'] == 'S', 'log_sigma'] = sigma4S
df = df0[df0['organism']== 'S']

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig2.suptitle('Syn '+str(vol)+'  Monocultures')

#format fig  
ax0.set_title('Syn'+str(vol)+ 'in 0 HOOH ') #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('Syn '+str(vol)+ 'in 400 HOOH ') #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)') #settign x axis label for graph 1
ax0.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)')#settign x axis label for graph 2 
ax1.set_ylabel('Cells(ml$^{-1}$)')  #setting y label for both subgraphs 
#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df0[df0['organism']== 'S']['time'],df0[df0['organism']== 'S']['avg1'],yerr=df0[df0['organism']== 'S']['std1'],color = c0, marker='o', label = 'avg1')
ax0.errorbar(df0[df0['organism']== 'S']['time'],df0[df0['organism']== 'S']['avg2'],yerr=df0[df0['organism']== 'S']['std2'], color = c0,marker='v', label = 'avg2')
ax0.errorbar(df0[df0['organism']== 'S']['time'],df0[df0['organism']== 'S']['abundance'],yerr=df0[df0['organism']== 'S']['sigma'], color = c0,marker='v', label = 'MEAN')

# graph 400 H assay even and odd avgs
ax1.errorbar(df4[df4['organism']== 'S']['time'],df4[df4['organism']== 'S']['abundance'],yerr=df4[df4['organism']== 'S']['std1'], color = c0,marker='o', label = 'avg1')
ax1.errorbar(df4[df4['organism']== 'S']['time'],df4[df4['organism']== 'S']['avg2'],yerr=df4[df4['organism']== 'S']['std2'],color = c0, marker='v', label = 'avg2')
ax1.errorbar(df4[df4['organism']== 'S']['time'],df4[df4['organism']== 'S']['abundance'],yerr=df4[df4['organism']== 'S']['sigma'],color = c0, marker='v', label = 'MEAN')

plt.legend()

#####################################################
#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 
#####################################################

#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
inits0 = pd.read_csv('../data/inits/syn_vol'+str(vol)+'_inits0.csv')

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
S0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+5})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
S0_mean = inits0['S0'][0]
N0_mean = inits0['N0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_0H,
                          parameter_names=['k1','k2','S0','N0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          S0 = S0_prior.copy(),
                          N0  = N0_prior.copy(),
                          t_steps=1000,
                          S = S0_mean,
                          N = N0_mean,
                            )
    return M

def mono_0H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2 = params[0], params[1]
    S,N = max(y[0],0),max(y[1],0),
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dSdt = (k2 * N /( (ksp) + N) )*S     
    dNdt =  - (k2 * N /( (ksp) + N) )*S
    return [dSdt,dNdt]

#find closest time 
def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)

#df0.loc[:,'log_abundance'] = np.log(10**df0.log_abundance)


# get_models
a0 = get_model(df) 

#broken here!!!!!!!!!!
# do fitting
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod0 = a0.integrate()

a0res = get_residuals(a0)  #is this using the best fit or just a first run???

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################

###### fig set up
fig3, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6)) #fig creationg of 1 by 2
fig3.suptitle('Syn '+str(vol)+ ' in 0 H Model') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.90, wspace = 0.45, left = 0.10, hspace = 0.20, bottom = 0.2)

ax0.semilogy()
ax0.set_title('Syn '+str(vol)+ '  dynamics ',fontsize = '16')
ax1.set_title('Model residuals',fontsize = '14')

ax0.set_xlabel('Time (days)',fontsize = '14')
ax0.set_ylabel('Cells(ml$^{-1}$)',fontsize = '14')
ax1.set_ylabel('Data S value',fontsize = '14')
ax1.set_xlabel('Residual',fontsize = '14')


l3 = ax0.legend(loc = 'lower right')
l3.draw_frame(False)


#graphing data from df to see 2 different biological reps represented

ax0.errorbar(df0[df0['organism']== 'S']['time'],df0[df0['organism']== 'S']['abundance'],yerr=df0[df0['organism']== 'S']['sigma'], color = c0,marker='v', label = 'MEAN S')
ax0.plot(mod0.time,mod0['S'],c='r',lw=1.5,label=' model best fit')
a0.plot_uncertainty(ax0,posteriors0,'S',1000)

ax1.scatter(a0res['res'], a0res['abundance'],label = '0H case')
#printing off graph

fig3.savefig('../figures/syn'+str(vol)+ '_0H_fits')


#########################################################
#graphing model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,3,figsize=[12,4])
#set titles and config graph 
#fig4.axes(fontsize = 10)


ax4[1].set_title(r'$S_0$', fontsize = 12)
ax4[2].set_title(r'$\mu_{max}$', fontsize = 12)
ax4[0].set_xlabel('Time (days)', fontsize = 12)
ax4[1].set_xlabel('Parameter Value', fontsize = 12)
ax4[1].set_ylabel('Frequency', fontsize = 12)
ax4[2].set_xlabel('Parameter Value', fontsize = 12)
ax4[2].set_ylabel('Frequency', fontsize = 12)
ax4[0].set_xlabel('Time (days)', fontsize = 14)
ax4[0].set_ylabel('Cells (ml$^{-1}$)', fontsize = 12)
ax4[0].tick_params(axis='x', labelsize=12)
ax4[0].tick_params(axis='y', labelsize=12)
ax4[1].tick_params(axis='x', labelsize=12)
ax4[1].tick_params(axis='y', labelsize=12)
ax4[2].tick_params(axis='x', labelsize=12)
ax4[2].tick_params(axis='y', labelsize=12)

#shift fig subplots
fig4.subplots_adjust(wspace = 0.3)

#graph data, model, and uncertainty 
ax4[0].plot(df0[df0['organism']== 'S']['time'], df0[df0['organism']== 'S']['abundance'],color = c0, marker='o',label = r'$Synechococcus$')
ax4[0].errorbar(df0[df0['organism']== 'S']['time'], df0[df0['organism']== 'S']['abundance'], yerr = df0[df0['organism']== 'S']['sigma'], color = c0, marker='o')
ax4[0].plot(mod0.time,mod0['S'],c='r',lw=1.5,label=' Model best fit')
a0.plot_uncertainty(ax4[0],posteriors0,'S',100)

# plot histograms of parameter search results 
ax4[1].hist(posteriors0.S0 , facecolor = c0)
ax4[2].hist(posteriors0.k2 , facecolor = c0)

ax4[0].semilogy()

for (ax,l) in zip(ax4,'abc'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#make legends
l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)

#show full graph 
fig4.savefig('../figures/syn '+str(vol)+ '_0H_odelib')


pframe = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe.to_csv('../data/inits/syn_vol'+str(vol)+ '_inits0.csv')

fig4.savefig('../figures/syn'+str(vol)+'_odelib0')


# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')


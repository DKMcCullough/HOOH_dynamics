
'''

name:   model_pro_batch1.py 

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
df_all = df_all.fillna(0)

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

df['log_sigma'] = 0.1

#slicing data into abiotic, biotic, and Pro only dataframes
df0 = df.loc[~ df['assay'].str.contains('4', case=False) & (df['Vol_number']== 1)]  #assay 0 H 
df4 = df.loc[(df['assay'].str.contains('4', case=False)) & (df['Vol_number']== 1)]


df = df0[df0['organism']== 'P']

c0 = 'g'

#####################################################
#plotting data and error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig2.suptitle('Pro  Monocultures')
fig2.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.20)


#format fig  
ax0.set_title('Pro in 0 HOOH ',fontsize = 16) #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('Pro in 400 HOOH ',fontsize = 13) #graph title for graph 2
ax1.semilogy()#setting y axis to be logged b/c cell data
ax0.set_xlabel('Time (days)',fontsize = 12) #settign x axis label for graph 1
ax0.set_ylabel('Cells (ml $^{-1}$)',fontsize = 12)  #setting y label for both subgraphs 
ax1.set_xlabel('Time (days)',fontsize = 12)#settign x axis label for graph 2 


#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df0[df0['organism']== 'P']['time'],df0[df0['organism']== 'P']['avg1'],yerr=df0[df0['organism']== 'P']['std1'], marker='o', label = 'Bio1')
ax0.errorbar(df0[df0['organism']== 'P']['time'],df0[df0['organism']== 'P']['avg2'],yerr=df0[df0['organism']== 'P']['std2'], marker='v', label = 'Bio2')
ax0.errorbar(df0[df0['organism']== 'P']['time'],df0[df0['organism']== 'P']['abundance'],yerr=df0[df0['organism']== 'P']['sigma'], marker='d', label = 'MEAN')

# graph 400 H assay even and odd avgs
ax1.errorbar(df4[df4['organism']== 'P']['time'],df4[df4['organism']== 'P']['avg1'],yerr=df4[df4['organism']== 'P']['std1'], marker='o', label = 'Bio1')
ax1.errorbar(df4[df4['organism']== 'P']['time'],df4[df4['organism']== 'P']['avg2'],yerr=df4[df4['organism']== 'P']['std2'], marker='v', label = 'Bio2')
ax1.errorbar(df4[df4['organism']== 'P']['time'],df4[df4['organism']== 'P']['abundance'],yerr=df4[df4['organism']== 'P']['sigma'], marker='d', label = 'MEAN')


l1 = ax0.legend(loc = 'lower right')
l1.draw_frame(False)

#fig2.savefig('../figures/pro1_data_0and4')

#####################################################
#plotting error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1,ax2)= plt.subplots(1,3,figsize = (11,5))
fig2.suptitle('Pro  in 0 H ')
fig2.subplots_adjust(right=0.95, wspace = 0.45, hspace = 0.20)


#format fig  
ax0.set_title('Pro Dyanimics ',fontsize = 15) #graph title for graph 1
ax1.set_title('Error Structure',fontsize = 13) #graph title for graph 2
ax2.set_title('LOG Error Structure',fontsize = 13) #graph title for graph 2

ax0.set_xlabel('Time (days)',fontsize = 12) #settign x axis label for graph 1
ax0.set_ylabel('Cells (ml $^{-1}$)',fontsize = 12)  #setting y label for both subgraphs 
ax0.semilogy()

ax1.set_xlabel('Abundance',fontsize = 12)#settign x axis label for graph 2 
ax1.set_ylabel('Variance',fontsize = 12)#settign x axis label for graph 2 
ax1.semilogx()
ax1.semilogy()

ax2.set_xlabel('Log Abundance',fontsize = 12)#settign x axis label for graph 2 
ax2.set_ylabel('Log Variance',fontsize = 12)#settign x axis label for graph 2 
ax2.semilogx()
ax2.semilogy()

#graph 0 H assay 
ax0.errorbar(df0[df0['organism']== 'P']['time'],df0[df0['organism']== 'P']['lavg1'],yerr=df0[df0['organism']== 'P']['stdlog1'], marker='o', label = 'Bio1')
ax0.errorbar(df0[df0['organism']== 'P']['time'],df0[df0['organism']== 'P']['lavg2'],yerr=df0[df0['organism']== 'P']['stdlog2'], marker='v', label = 'Bio2')
ax0.errorbar(df0[df0['organism']== 'P']['time'],df0[df0['organism']== 'P']['log_abundance'],yerr=df0[df0['organism']== 'P']['log_sigma'], marker='d', label = 'MEAN')

# raw data vs sigma
ax1.scatter(df0[df0['organism']== 'P']['avg1'],df0[df0['organism']== 'P']['std1'], marker='o', label = 'Bio1')
ax1.scatter(df0[df0['organism']== 'P']['avg2'],df0[df0['organism']== 'P']['std2'], marker='v', label = 'Bio2')
ax1.scatter(df0[df0['organism']== 'P']['abundance'],df0[df0['organism']== 'P']['sigma'], marker='d', label = 'MEAN')

# log data vs sigma
ax2.scatter(df0[df0['organism']== 'P']['lavg1'],df0[df0['organism']== 'P']['stdlog1'], marker='o', label = 'Log_Bio1')
ax2.scatter(df0[df0['organism']== 'P']['lavg2'],df0[df0['organism']== 'P']['stdlog2'], marker='v', label = 'Log_Bio2')
ax2.scatter(df0[df0['organism']== 'P']['log_abundance'],df0[df0['organism']== 'P']['log_sigma'], marker='d', label = 'Log_MEAN')

l1 = ax0.legend(loc = 'lower right')
l1.draw_frame(False)

#fig2.savefig('../figures/pro1_0_data_trends')

#####################################################
#plotting error within biological reps 
#####################################################
# fig set up and main title 
fig2, (ax0,ax1,ax2)= plt.subplots(1,3,figsize = (11,5))
fig2.suptitle('Pro  in 400 H ')
fig2.subplots_adjust(right=0.95, wspace = 0.45, hspace = 0.20)

#format fig  
ax0.set_title('Pro Dyanimics ',fontsize = 15) #graph title for graph 1
ax0.semilogy() #setting y axis to be logged b/c cell data
ax1.set_title('Error Structure',fontsize = 13) #graph title for graph 2
ax2.set_title('LOG Error Structure',fontsize = 13) #graph title for graph 2


ax0.set_xlabel('Time (days)',fontsize = 12) #settign x axis label for graph 1
ax0.set_ylabel('Cells (ml $^{-1}$)',fontsize = 12)  #setting y label for both subgraphs 
ax0.semilogy()

ax1.set_xlabel('Abundance',fontsize = 12)#settign x axis label for graph 2 
ax1.set_ylabel('Variance',fontsize = 12)#settign x axis label for graph 2 
ax1.semilogx()
ax1.semilogy()

ax2.set_xlabel('Log Abundance',fontsize = 12)#settign x axis label for graph 2 
ax2.set_ylabel('Log Variance',fontsize = 12)#settign x axis label for graph 2 
ax2.semilogx()
ax2.semilogy()
#graph dataframe of even or odd avgs (for tech reps) to give avg of total bioreps 

#graph 0 H assay even and odd avgs 
ax0.errorbar(df4[df4['organism']== 'P']['time'],df4[df4['organism']== 'P']['avg1'],yerr=df4[df4['organism']== 'P']['std1'], marker='o', label = 'Bio1')
ax0.errorbar(df4[df4['organism']== 'P']['time'],df4[df4['organism']== 'P']['avg2'],yerr=df4[df4['organism']== 'P']['std2'], marker='v', label = 'Bio2')
ax0.errorbar(df4[df4['organism']== 'P']['time'],df4[df4['organism']== 'P']['abundance'],yerr=df4[df4['organism']== 'P']['sigma'], marker='d', label = 'MEAN')

# graph 400 H assay even and odd avgs
ax1.scatter(df4[df4['organism']== 'P']['avg1'],df4[df4['organism']== 'P']['std1'], marker='o', label = 'Bio1')
ax1.scatter(df4[df4['organism']== 'P']['avg2'],df4[df4['organism']== 'P']['std2'], marker='v', label = 'Bio2')
ax1.scatter(df4[df4['organism']== 'P']['abundance'],df4[df4['organism']== 'P']['sigma'], marker='d', label = 'MEAN')

ax2.scatter(df4[df4['organism']== 'P']['lavg1'],df4[df4['organism']== 'P']['stdlog1'], marker='o', label = 'Log_Bio1')
ax2.scatter(df4[df4['organism']== 'P']['lavg2'],df4[df4['organism']== 'P']['stdlog2'], marker='v', label = 'Log_Bio2')
ax2.scatter(df4[df4['organism']== 'P']['log_abundance'],df4[df4['organism']== 'P']['log_sigma'], marker='d', label = 'Log_MEAN')


l1 = ax0.legend(loc = 'upper right')
l1.draw_frame(False)

#fig2.savefig('../figures/pro1_4_data_trends')

#####################################################
#####################################################

#   model param and state variable set up 
# modeling abiotic HOOH via SH and deltaH and H0 

#####################################################
#####################################################
#reading in csv file with inititl guesses for all parameter values ( SH, deltah, H0)
inits0 = pd.read_csv("../data/inits/pro_MIT9215_inits0_1.csv")

#setting how many MCMC chains you will run 
nits = 10000 # number of iterations - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
snames = ['P','N'] #order must match all further model mentions (same fro params) 

#sigma for param search
pw = 1   


#setting param prior guesses and inititaing as an odelib param class in odelib
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':2000})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
#setting state variiable  prior guess
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2e+8})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#set mean for model via inits
P0_mean = inits0['P0'][0]
N0_mean = inits0['N0'][0]

#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_0H,
                          parameter_names=['k1','k2','P0','N0'],
                          state_names = snames,
                          dataframe=df,
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          P0 = P0_prior.copy(),
                          N0  = N0_prior.copy(),
                          t_steps=1000,
                          P = P0_mean,
                          N = N0_mean,
                            )
    return M

def mono_0H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2 = params[0], params[1]
    #P,N = max(y[0],0),max(y[1],0),
    P,N = y[0],y[1]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P     
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    return [dPdt,dNdt]

#find closest time 
def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)


# get_model of df using function
a0 = get_model(df) 

# do fitting of model witing MCMC method
posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod0 = a0.integrate()

#get residuals between best model fit and data 
a0res = get_residuals(a0)  

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################

###### fig set up
fig3, (ax0,ax1)= plt.subplots(1,2,figsize = (9,4)) #fig creationg of 1 by 2
#fig3.suptitle('Pro in 0 H Model',fontsize = '16') #setting main title of fig

####### fig config and naming 

fig3.subplots_adjust(right=0.9, wspace = 0.45, hspace = 0.20, left = 0.1)

ax0.semilogy()
ax0.set_title('Pro dynamics ',fontsize = '14')
ax1.set_title('Model residuals',fontsize = '14')

ax0.set_xlabel('Time (days)',fontsize = '14')
ax0.set_ylabel('Cells(ml$^{-1}$)',fontsize = '14')
ax1.set_ylabel('Data P value',fontsize = '14')
ax1.set_xlabel('Residual',fontsize = '14')



#model and residuals
ax0.errorbar(df[df['organism']== 'P']['time'],df[df['organism']== 'P']['abundance'],yerr=df[df['organism']== 'P']['sigma'],c=c0, marker='d', label = 'Pro Data MEAN')
ax0.plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' model best fit')
a0.plot_uncertainty(ax0,posteriors0,'P',100)

ax1.scatter(a0res['res'], a0res['abundance'],c=c0, label = '0H case')

#printing off graph
l3 = ax0.legend(loc = 'lower right')
l3.draw_frame(False)

#fig3.savefig('../figures/pro1_odelib0_fit')

#########################################################
#graphing model vs data and params histograms 
#########################################################

# set up graph
fig4,ax4 = plt.subplots(1,3,figsize=[12,4])
#set titles and config graph 
fig4.subplots_adjust(wspace = 0.40, bottom=0.2)
ax4[0].set_ylim([2e+5, 4e+6])

#ax4[0].set_title('Pro  dynamics', fontsize = 16)
#plt.text(0.5, 1.08, 'Pro  dynamics',horizontalalignment='left',fontsize=14, transform = ax4.transAxes)
ax4[0].set_xlabel('Time (days)', fontsize = 12)
ax4[1].set_xlabel('Initial cell density \n ($P_{i,0}$, cells mL$^{-1}$)', fontsize = 12)
ax4[1].set_ylabel('Probability density (x10$^{-6}$)', fontsize = 12)
ax4[2].set_xlabel('Growth rate ($\mu_i$, day$^{-1}$)', fontsize = 12)
ax4[2].set_ylabel('Probability density', fontsize = 12)
ax4[0].set_ylabel('Cells (mL$^{-1}$)', fontsize = 12)
ax4[0].tick_params(axis='x', labelsize=12)
ax4[0].tick_params(axis='y', labelsize=12)
ax4[1].tick_params(axis='x', labelsize=12)
ax4[1].tick_params(axis='y', labelsize=12)
ax4[2].tick_params(axis='x', labelsize=12)
ax4[2].tick_params(axis='y', labelsize=12)

#shift fig subplots
fig4.subplots_adjust(wspace = 0.3)
ax4[0].semilogy()

for (ax,l) in zip(ax4,'abc'):
    ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)

#graph data, model, and uncertainty 
ax4[0].errorbar(df[df['organism']== 'P']['time'],df[df['organism']== 'P']['abundance'],yerr=df[df['organism']== 'P']['sigma'], c=c0, marker='d', label = r'$Prochlorococcus$')
ax4[0].plot(mod0.time,mod0['P'],c='r',lw=1.5,label='Model best fit')
a0.plot_uncertainty(ax4[0],posteriors0,'P',100)

# plot histograms of parameter search results 
ax4[1].hist(posteriors0.P0,density=True, color = c0)
ax4[2].hist(posteriors0.k2,density=True, color = c0)

#format legend 
l4 = ax4[0].legend(loc = 'lower right', fontsize = 9)
l4.draw_frame(False)
#show full graph 

# rescale
ticks = ax4[1].get_yticks()
ax4[1].set_yticklabels([f"{tick * 1e+6:.0f}" for tick in ticks])

fig4.savefig('../figures/pro1_odelib0_params',bbox_inches='tight')

#TRACE plot for growth params
fig5,ax5 = plt.subplots(1,2,sharex=True,figsize=[8,4]) #make plot
fig5.suptitle('Trace plots for Params ', fontsize = 14) #set main title 
fig5.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8, wspace=0.45, hspace=0.2) #shift white space for better fig view

ax5[0].set_title('\u03B1', fontsize = 14)
ax5[1].set_title('\u03BC', fontsize = 14)

ax5[0].set_ylabel('\u03B1 value', fontsize = 12)
ax5[0].set_xlabel('Model iteration', fontsize = 12)
ax5[1].set_ylabel('\u03BC value', fontsize = 12)
ax5[1].set_xlabel('Model iteration', fontsize = 12)

#ax3[:,:].set_yscale('log')


#graphing iteration number vs parameter numbert logged 
ax5[0].scatter(posteriors0.iteration,posteriors0.k1,color = c0)
ax5[1].scatter(posteriors0.iteration,posteriors0.k2,color = c0)


#print out plot
#fig3.savefig('../figures/Pro1_0_TRACE')



#update inits file withg best model params
pframe = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
pframe.to_csv('../data/inits/pro_MIT9215_inits0_1.csv')



# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')

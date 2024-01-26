'''
HEPES_3_way_solved.py

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Using Morris_et_al_2011 data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

import scipy
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   

import ODElib
import random as rd
import sys



###############################

#  Importing the data frame

###############################

 #UH18301 Pro in 3.75 mM HEPES or Taps buffer. High light (24uC in a Sunbox -  noon maximum of about 250 quanta m) 

df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'Buffers_Morris_2011_f1c,d', header  = 0)

   #'renaming column to make it callable by 'times'
df_all = df_all.fillna(0)

#df_all = pd.read_csv('../data/Buffers_Morris_2011_f1.csv')

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time (days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#creating log and stats for data 
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['sigma'] =  np.std(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)


##############################

#   data arrays

###############################


orgs = df_all['organism'].unique()

assays = df_all['ID'].unique()
nassays = assays.shape[0]
treats  = df_all['Treatment'].unique()
ntreats = treats.shape[0]

#HOOHs_nm = 1000**np.array(HOOH_data)   #TRIED TO GET IN nM to match better with other things 
#print(HOOH_data)
colors = ('r','g')
markers = ('o','*')

##############################

#    graphing the data 

##############################


for a,n in zip(assays,range(nassays)):
    fig1,(ax1)= plt.subplots(2,1, figsize = (12,8))
    for t,nt in zip(treats,range(ntreats)): 
        df = df_all[(df_all['ID'] == a)]
        count = nt
        df = df[(df['ID']==a) & (df['Treatment']==t)]
        pdf = df[df['organism'] == 'P']
        hdf = df[df['organism'] == 'H']
        ax1[0].errorbar(pdf['time'], pdf['abundance'], yerr = pdf['sigma'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        ax1[1].errorbar(hdf['time'], hdf['abundance'], yerr = hdf['sigma'],marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        fig1.suptitle('Dynamics of '+ str(a))
        ax1[1].set_ylabel('HOOH concentration (\u03BCM)')
        ax1[0].set_ylabel('Pro abundance (cell/ml)')
        fig1.supxlabel('Time (days)')
        l1 = ax1[0].legend(loc = 'lower right', prop={"size":14}) 
        l1.draw_frame(False)#print(df)
        ax1[0].semilogy()
        ax1[1].semilogy()
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)

    plt.show()
    fig1.savefig('../figures/Hproduction_'+str(a)+'_data.png')


#############################
#split into just hepes df for checking different types od modeling directly against each other
##############################

df = df_all[(df_all['ID'] =='UH18301') & (df_all['Treatment']=='HEPES') & (df_all['organism'] == 'H')]



fig2,(ax2) = plt.subplots(figsize = (7,6))
fig2.suptitle('HOOH production from HEPES buffer', size = 16)
ax2.set_ylabel('HOOH concentration (\u03BCM)', size = 14)

ax2.set_xlabel('Time (days)', size = 14)

ax2.errorbar(x = df['time'],y = df['abundance'],yerr= df['sigma'], c='purple', marker = 'o', linestyle = ':',label='HEPES mean data')


####################################

#analytical solution

####################################

#initial values and creating time array

delta = 0.2
S_HOOH = 690
step = 0.05 #delta t
ndays = 7
H0 = 40
times = np.linspace(0,ndays,int(ndays/step))

def f(t, S_HOOH, delta):
    H = (S_HOOH/delta)*(1-(np.exp(-delta*t))) + (H0*(np.exp(-delta*t)))
    return H

#calculating Hs from analutical solution over a range of ts 
Hs = f(times,S_HOOH,delta)

#plot Hs on fig with data
ax2.plot(times,Hs,c='c',linestyle = '-',label='Model via Analytical Solution')

#############################

#Euler's Integration

############################

HsEuler = np.array([]) 
H = 40
t0 = 0

for t in times: 
	HsEuler = np.append(HsEuler,H)
	dHdt = S_HOOH - delta*H
	H = H + dHdt*step
	
    #plot Eulers soluton of Hs
ax2.plot(times,HsEuler,c='g',linestyle = '--',label = "Model via Euler's Aproximation")#,label = "Euler's Aproximation")


####################################

#ODE int

####################################


from scipy.integrate import odeint

def HsODEint(H,t):
	dHdt = S_HOOH-delta*H
	#print(dHdt,t,H)
	return dHdt


ode_solutions = odeint(HsODEint,3.5,times)


plt.plot(times,ode_solutions,c='dodgerblue', linestyle = ':', label = 'Model via Odeint Approximation')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)



#####################################################

#ODElib

#####################################################

inits4 = pd.read_csv("../data/inits/hepes_3way.csv")


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



#find closest time to that of data in best model and comparing values
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
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':300})

priors = {'deltah' :deltah_prior,'Sh' : Sh_prior,'H0' : H0_prior} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = inits4['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 10000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a4 = get_model(df) 

#broken here!!!!!!!!!!
# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()
#calculating residuals of model
a4res = get_residuals(a4)

#plot model and search
plt.plot(mod4.time,mod4['H'],c='r',lw=1.5,label = 'Model best fit via Odelib') #best model fit of 0 H assay
a4.plot_uncertainty(plt,posteriors4,'H',100)



#set full figure legend
plt.legend(loc = 'lower right')
plt.show()

	
fig2.savefig('../figures/HEPES_3way')


pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/hepes_3way.csv')


'''

################
############
######
#taps 
######
############
################

df = df_all[(df_all['ID'] =='UH18301') & (df_all['Treatment']=='TAPS') & (df_all['organism'] == 'H')]



fig2,(ax2) = plt.subplots(figsize = (7,6))
fig2.suptitle('HOOH production from TAPS buffer', size = 16)
ax2.set_ylabel('HOOH concentration (\u03BCM)', size = 14)

ax2.set_xlabel('Time (days)', size = 14)

ax2.errorbar(x = df['time'],y = df['abundance'],yerr= df['sigma'], c='purple', marker = 'o', linestyle = ':',label='TAPS mean data')


####################################

#analytical solution

####################################

#initial values and creating time array


delta = 0.2
S_HOOH = 32
step = 0.05 #delta t
ndays = 7
H0 = 25
times = np.linspace(0,ndays,int(ndays/step))

def f(t, S_HOOH, delta):
    H = (S_HOOH/delta)*(1-(np.exp(-delta*t))) + (H0*(np.exp(-delta*t)))
    return H

#calculating Hs from analutical solution over a range of ts 
Hs = f(times,S_HOOH,delta)

#plot Hs on fig with data
ax2.plot(times,Hs,c='c',linestyle = '-',label='Model via Analytical Solution')

#############################

#Euler's Integration

############################

HsEuler = np.array([]) 
H = 1
t0 = 0

for t in times: 
    
	#HsEuler = np.append(HsEuler,H)
    if t == 0: 
        H = 25
    else: 
        dHdt = S_HOOH - delta*H 
        H =  H +  dHdt*step 
    HsEuler = np.append(HsEuler,H)


ax2.plot(times,HsEuler,c='g',linestyle = '--',label = "Model via Euler's Aproximation")#,label = "Euler's Aproximation")


####################################

#ODE int

####################################


from scipy.integrate import odeint

def HsODEint(H,t):
	dHdt = S_HOOH-delta*H 
	#print(dHdt,t,H)
	return dHdt


ode_solutions = odeint(HsODEint,H0,times)


plt.plot(times,ode_solutions,c='dodgerblue', linestyle = ':', label = 'Model via Odeint Approximation')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)



#####################################################

#ODElib

#####################################################

inits4 = pd.read_csv("../data/inits/taps_3way.csv")


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



#find closest time to that of data in best model and comparing values
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
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':300})

priors = {'deltah' :deltah_prior,'Sh' : Sh_prior,'H0' : H0_prior} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = inits4['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 10000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a4 = get_model(df) 

#broken here!!!!!!!!!!
# do fitting
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
mod4 = a4.integrate()
#calculating residuals of model
a4res = get_residuals(a4)

#plot model and search
plt.plot(mod4.time,mod4['H'],c='r',lw=1.5,label = 'Model best fit via Odelib') #best model fit of 0 H assay
a4.plot_uncertainty(plt,posteriors4,'H',100)



#set full figure legend
plt.legend(loc = 'lower right')
plt.show()

	
fig2.savefig('../figures/TAPS_3way')


pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
pframe.to_csv('../data/inits/taps_3way.csv')

'''



print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')

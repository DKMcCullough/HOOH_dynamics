
'''

name:   model_pro_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model and graph 0 H Pro assay and model of said biomass via odelib

working on: ln of data in df for uncertainty, loop for 0 and 400 using different init files? (need H connected to 0 H first) 

'''









import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ODElib
import random as rd
import sys


plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'


#########
#reading in data
############
df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#slicing data

df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 

#setting working df
df = df_P

#making avg columns of technical reps (std hereo nly for graphing, not logged here)
df['avg1'] = df[['rep1', 'rep2']].mean(axis=1)
df['avg2'] = df[['rep3', 'rep4']].mean(axis=1)
df['std1'] = df[['rep1', 'rep2']].std(axis=1)
df['std2'] = df[['rep3', 'rep4']].std(axis=1)


df.rename(columns = {'avg1':'abundance'}, inplace = True)


df0 = df.loc[~ df['assay'].str.contains('4', case=False)] 
df4 = df.loc[df['assay'].str.contains('4', case=False)] 


fig2, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig2.suptitle('Pro  Monocultures')

ax0.errorbar(df0['time'],df0['abundance'],yerr=df0['std1'], marker='o', label = 'avg1')
ax0.errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='v', label = 'avg2')
ax0.set_title('Pro in 0 HOOH ')
ax0.semilogy()
ax1.errorbar(df4['time'],df4['abundance'],yerr=df4['std1'], marker='o', label = 'avg1')
ax1.errorbar(df4['time'],df4['avg2'],yerr=df4['std2'], marker='v', label = 'avg2')
ax1.set_title('Pro in 400 HOOH ')
ax1.semilogy()
ax0.set_xlabel('Time (days)')
ax0.set_ylabel('Cells(ml$^{-1}$)')
ax1.set_xlabel('Time (days)')


#########
# modeling abiotic for SH and deltaH and H0 info
#########


inits0 = pd.read_csv("../data/inits/pro9215_inits0.csv")


# state variable names
snames = ['P','N','H']

# define priors
pw = 1

Qnp = int((9.4e-15*(1/(14.0))*1e+9))  #Nitrogen Quota for Pro from Bertilison 

Qnp_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.0000020})
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.0000002})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
dp_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.002})
rho_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
SN_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':4})
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1})
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+6})
N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})



nits = 100 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS

P0_mean = 100000
N0_mean = 900000

H0_mean = 80



def set_best_params(model,posteriors,snames):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    bestchain = posteriors.iloc[im]["chain#"]
    posteriors = posteriors[posteriors["chain#"]==bestchain]
    model.set_parameters(**posteriors.loc[im][a1.get_pnames()].to_dict())
    model.set_inits(**{o:posteriors.loc[im][a1.get_pnames()].to_dict()[o+'0'] for o in ['P']})
#    model.set_inits(**{o:posteriors.loc[im][a1.get_pnames()].to_dict()[o+'0'] for o in ['H']})

def plot_uncertainty(ax,model,posteriors,ntimes):
    for a in range(ntimes):
        im = rd.choice(posteriors.index)
        model.set_inits(**{'P':posteriors.loc[im][model.get_pnames()].to_dict()['P0']})
        #model.set_inits(**{'H':posteriors.loc[im][model.get_pnames()].to_dict()['H0']})
        model.set_parameters(**posteriors.loc[im][model.get_pnames()].to_dict())
        mod = model.integrate()
        ax.plot(mod.time,mod['P'],c=str(0.8),lw=1,zorder=1)
        #ax.plot(mod.time,mod['H'],c=str(0.8),lw=1,zorder=1)

#Ksp or k1?!?

def get_model(df):
    a1=ODElib.ModelFramework(ODE=mono_0H,
                          parameter_names=['deltah','Sh','rho','Qnp','SN','k1','k2','dp','P0','N0','H0'],
                          state_names = snames,
                          dataframe=df,
                          deltah = deltah_prior.copy(),
                          Sh = Sh_prior.copy(),
                          rho = rho_prior.copy(),
                          Qnp = Qnp_prior.copy(),
                          SN = SN_prior.copy(),
                          k1 = k1_prior.copy(),
                          k2 = k2_prior.copy(),
                          dp = dp_prior.copy(),
                          P0 = P0_prior.copy(),
                          N0  = N0_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          P = P0_mean,
                          N = N0_mean,
                          H = H0_mean
                            )
    return a1




#get k2, ksp, dp fit here and maybe rho and N0 or SN too?
def mono_0H(y,t,params): #no kdam or phi here (or make 0)
    deltah,Sh,rho,Qnp,SN,k1,k2,dp = params[0], params[1], params[2], params[3], params[4], params[5],params[6],params[7]
    P,N,H = y[0],y[1],y[2]
    ksp=k2/k1
    #print(P)
    dPdt = (k2 * N /( (ksp) + N) )*P - (dp *P)     
    dNdt =  SN - ((k2 * N /( (ksp) + N) )*P* Qnp) - rho*N    
    dHdt = Sh - deltah*H  #phi being P cell-specific detox rate
    return [dPdt,dNdt,dHdt]


df0.loc[:,'log_abundance'] = np.log(10**df0.log_abundance)

# get_models
a1 = get_model(df0) 


# do fitting
posteriors1 = a1.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1,static_parameters=set(['Qnp']))
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# set best params
set_best_params(a1,posteriors1,snames)

# run model with optimal params
mod0 = a1.integrate()



####################
# graphing model
#####################


fig3, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6))
fig3.suptitle('Pro in 0 H Model')

ax0.errorbar(df0['time'],df0['abundance'],yerr=df0['std1'], marker='o', label = 'avg1')
ax0.errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='o', label = 'avg2')

ax0.plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax0,a1,posteriors1,100)
ax0.semilogy()
ax0.set_title('Pro dynamics ')
l3 = ax0.legend(loc = 'lower right')
l3.draw_frame(False)

ax1.scatter(posteriors1.iteration, posteriors1.chi)
ax1.set_title('Model error ')

fig3.subplots_adjust(right=0.85, wspace = 0.25, hspace = 0.30)

ax0.set_xlabel('days')
ax0.set_ylabel('cell concentration')
ax1.set_ylabel('chi')
ax1.set_xlabel('iteration number ')

plt.show()


#########################################################
#graphing model and params 
#########################################################
'''
# pro model graph
fig4,ax4 = plt.subplots(1,7,figsize=[20,7])
fig4.suptitle('Monoculture parameters in 0 HOOH ')
ax4[0].plot(df0.time,df0.abundance, marker='o',label = 'Pro Mono - 0 H ')
ax4[0].plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax4[0],a1,posteriors1,100)

l4 = ax4[0].legend(loc = 'upper left')
l4.draw_frame(False)

# plot histograms
ax4[1].hist(posteriors1.dp)
ax4[2].hist(posteriors1.k1)
ax4[3].hist(posteriors1.k2)
ax4[4].hist(posteriors1.rho)
ax4[5].hist(posteriors1.Sh)
ax4[6].hist(posteriors1.deltah)


ax4[1].set_title('dp')
ax4[2].set_title('k1')
ax4[3].set_title('k2')
ax4[4].set_title('rho')
ax4[5].set_title('Sh')
ax4[6].set_title('deltah')


ax4[3].set_xlabel('Frequency')


fig4.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)

plt.show()


fig4.savefig('../figures/pro_odelib0')
'''

print("I'm done bro")




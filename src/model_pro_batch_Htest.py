import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ODElib
import random as rd

df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df = df_all

df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  
df_co = df_all.loc[df_all['assay'].str.contains('coculture', case=False)].copy()  
df_mono = df_all.loc[~df_all['assay'].str.contains('coculture', case=False)].copy()  

#df = df_mono 

#df['log_abundance'] = np.nanmean(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))
#df['log_stdv'] = np.std(np.log(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4']]],axis=0))



df_P = df_mono.loc[df_mono['organism'].str.contains('P', case=False)].copy() 
df_S = df_mono.loc[df_mono['organism'].str.contains('S', case=False)].copy() 


df = df_mono

plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'

strains = df_mono['strain'].unique()
nstrains = strains.shape[0]
colors = ('orange', 'r', 'green', 'c', 'purple', 'k')
#need to slice by Vol number !!! (2 cat + Syns)
fig1, (ax0,ax1)= plt.subplots(1,2,figsize = (16,9))
fig1.suptitle('Pro and Syn Monocultures')

for (S,si) in zip(strains,range(nstrains)):   
    count = si
    df = df_mono[((df_mono['strain']==S))].copy()
    df['log1'] = np.log(df['rep1'])
    df['log2'] = np.log(df['rep2'])
    df['log3'] = np.log(df['rep3'])
    df['log4'] = np.log(df['rep4'])
    df['avg1'] = df[['log1', 'log3']].mean(axis=1)
    df['avg2'] = df[['log2', 'log4']].mean(axis=1)
    df['std1'] = df[['log1', 'log3']].std(axis=1)
    df['std2'] = df[['log2', 'log4']].std(axis=1)
    df0 = df[((df['assay']=='plus_0'))].copy()
    df400 = df[((df['assay']=='plus_400'))].copy()
    ax0.errorbar(df0['time'],df0['avg1'],yerr=df0['std1'], marker='o',color = colors[count], label = str(S)+' avg 1')
    ax0.errorbar(df0['time'],df0['avg2'],yerr=df0['std2'], marker='v',color = colors[count], label = str(S)+' avg 2')
    ax1.errorbar(df400['time'],df400['avg1'], yerr=df400['std1'],marker='o', color = colors[count], label = str(S)+' avg 1')
    ax1.errorbar(df400['time'],df400['avg2'], yerr=df400['std2'],marker='v', color = colors[count], label = str(S)+' avg 2')
    ax1.set_ylim(6.5,15)
    ax0.set_ylim(6.5,15)
    l1  = ax0.legend(loc = 'lower left')
    l1.draw_frame(False)
    
# make space on the right for annotation (e.g. ROS=0, etc.)
fig1.subplots_adjust(right=0.85, wspace = 0.20, hspace = 0.25)

# titles and labels 
ax0.set_title('Monocultures in 0 HOOH')
ax1.set_title('Monocultures in 400 HOOH')

ax0.set_xlabel('Time (days)')
ax1.set_xlabel('Time (days)')
ax0.set_ylabel('log(Cells(ml$^{-1}$))')
ax1.set_ylabel('log(Cells(ml$^{-1}$))')

#fig1.savefig('../figures/mono_all_graphed')



#########
# modeling abiotic for SH and deltaH and H0 info
#########

df = df_P
#df['log_abundance'] = np.log(df['raw_abundance'])
#df['log_sigma'] = np.std(df['HOOH_stdv'])

df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['log4'] = np.log(df['rep4'])
df['avg1'] = df[['log1', 'log3']].mean(axis=1)
df['avg2'] = df[['log2', 'log4']].mean(axis=1)
df['std1'] = df[['log1', 'log3']].std(axis=1)
df['std2'] = df[['log2', 'log4']].std(axis=1)

df = df.rename(columns={"avg1": "abundance"})

df0 = df.loc[~ df['assay'].str.contains('4', case=False)] 
df4 = df.loc[df['assay'].str.contains('4', case=False)] 

inits0 = pd.read_csv("../data/inits/pro9215_inits0.csv")
inits0 = pd.read_csv("../data/inits/pro9215_inits4.csv")

# state variable names
snames = ['H','P','N']

# define priors
pw = 1
Qnp = int((9.4e-15*(1/(14.0))*1e+9))  #Nitrogen Quota for Pro from Bertilison 

Qnp_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
dp_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
rho_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
SN_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':6})
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':6})


N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})
P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})



nits = 10 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS

N0_mean = 500
P0_mean = df.loc[df['time'] == 0, 'abundance'].iloc[0]
H0_mean = 400



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
                          parameter_names=['deltah','Sh','rho','Qnp','SN','k1','k2','dp','H0','N0','P0'],
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
                          H0  = H0_prior.copy(),
                          N0  = N0_prior.copy(),
                          P0 = P0_prior.copy(),
                          t_steps=10000,
                          N = N0_mean,
                          P = P0_mean,
                          H = H0_mean,
                            )
    return a1




#get k2, ksp, dp fit here and maybe rho and N0 or SN too?
def mono_0H(y,t,params): #no kdam or phi here (or make 0)
    deltah,Sh,rho,Qnp,SN,k1,k2,dp = params[0], params[1], params[2], params[3], params[4], params[5],params[6],params[7]
    P,N,H = y[0],y[1],y[2]
    ksp=k2/k1
    dPdt = (k2 * N /( (ksp) + N) )*P - (dp *P)     
    dNdt =  SN - ((k2 * N /( (ksp) + N) )*P* Qnp) - rho*N    
    dHdt = Sh - deltah*H  #phi being P cell-specific detox rate
    return [dPdt,dNdt,dHdt]



#params = (deltah,Sh,rho,Qnp,SN,ksp,k2,dp,phip,kdam)


'''

def mono_400H(y,t,params): #kdam and phi are being fit, all else should come from mono_0H fits
    deltah,Sh,rho,Qnp,SN,ksp,k2,dp,phip,kdam = params[0], params[1], params[2], params[3],params[4], params[5], params[6], params[7], params[8],params[9]
    P,N,H = y[0],y[1],y[2]
    dPdt = (k2 * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P     
    dNdt =  SN - ((k2 * N /( (ksp) + N) )*P* Qnp) - rho*N    
    dHdt = Sh - deltah*H -phip*H*P #phi being S cell-specific detox rate
    return [dPdt,dNdt,dHdt]

'''


# get_models
a1 = get_model(df0) 
#a4 = get_model(df4) 
 


# do fitting
posteriors1 = a1.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1)
#posteriors4 = a4.MCMC(chain_inits=chain_inits,iterations_per_chain=nits,cpu_cores=1)

# set best params
set_best_params(a1,posteriors1,snames)
#set_best_params(a4,posteriors4,snames)
# run model with optimal params
mod0 = a1.integrate()

#mod4 = a4.integrate()




#########################################################
# modeling
#########################################################

# pro model graph
fig4,ax4 = plt.subplots(2,3,figsize=[8,4.5])
ax4[0,0].plot(df0.time,df0.abundance, marker='o',label = 'Pro Mono - 0 H ')
ax4[0,0].plot(mod0.time,mod0['P'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax4[0,0],a1,posteriors1,100)
ax4[1,0].plot(df4.time,df4.abundance, marker='o',label = 'Pro Mono - 400 H ')
#ax4[1,0].plot(mod4.time,mod4['H'],c='r',lw=1.5,label=' model best fit')
#plot_uncertainty(ax4[1,0],a4,posteriors1,100)

ax4[0,0].legend(loc = 'upper left')
ax4[1,0].legend(loc = 'upper left')



# plot histograms
ax4[0,1].hist(posteriors1.dp)
ax4[0,2].hist(posteriors1.k1)
#ax4[1,1].hist(posteriors4.Sh)
#ax4[1,2].hist(posteriors4.deltah)


plt.show()


fig4.savefig('../figures/pro_odelib')




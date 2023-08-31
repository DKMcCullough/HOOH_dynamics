import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ODElib
import random as rd

#####################################################
# read in data and formatting
#####################################################

df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df = df_all

df_abiotic = df_all.loc[df_all['assay'].str.contains('abiotic', case=False)].copy()  

#########
# modeling abiotic for SH and deltaH and H0 info
#########

df = df_abiotic
df['log_abundance'] = np.log(df['raw_abundance'])
#df['log_sigma'] = np.std(df['HOOH_stdv'])
df['log_sigma'] = 0.1 # made up number
df = df.rename(columns={"raw_abundance": "abundance"})

df0 = df.loc[~ df['assay'].str.contains('4', case=False)] 
df4 = df.loc[df['assay'].str.contains('4', case=False)] 

## Reading in inits for 0 and 400 models
inits0 = pd.read_csv("../data/inits/abiotic0.csv")
inits4 = pd.read_csv("../data/inits/abiotic4.csv")

def set_best_params(model,posteriors,snames):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    bestchain = posteriors.iloc[im]["chain#"]
    posteriors = posteriors[posteriors["chain#"]==bestchain]
    model.set_parameters(**posteriors.loc[im][a1.get_pnames()].to_dict())
    model.set_inits(**{o:posteriors.loc[im][a1.get_pnames()].to_dict()[o+'0'] for o in ['H']})

def plot_uncertainty(ax,model,posteriors,ntimes):
    for a in range(ntimes):
        im = rd.choice(posteriors.index)
        model.set_inits(**{'H':posteriors.loc[im][model.get_pnames()].to_dict()['H0']})
        model.set_parameters(**posteriors.loc[im][model.get_pnames()].to_dict())
        mod = model.integrate()
        ax.plot(mod.time,mod['H'],c=str(0.8),lw=1,zorder=1)


def get_model(df):
    a1=ODElib.ModelFramework(ODE=abiotic,
                          parameter_names=['deltah','Sh', 'H0'],
                          state_names = snames,
                          dataframe=df,
                          deltah = deltah_prior.copy(),
                          Sh = Sh_prior.copy(),
                          H0  = H0_prior.copy(),
                          t_steps=1000,
                          H = H0_mean,
                         )
    return a1


def abiotic(y,t,params):
    deltah,Sh = params[0], params[1]
    H = y[0]
    dHdt = Sh - deltah*H 
    return [dHdt]


pw = 1

deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})

Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':6})

H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1e+5})


H0_mean = df.loc[df['time'] == 0, 'abundance'].iloc[0]




# state variable names
snames = ['H']

# get_models
a1 = get_model(df0) 
a4 = get_model(df4) 
 



# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 100000

# do fitting
posteriors1 = a1.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1)
posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1)


# set best params
set_best_params(a1,posteriors1,snames)
set_best_params(a4,posteriors4,snames)

# run model with optimal params
mod0 = a1.integrate()
mod4 = a4.integrate()



#########################################################
# modeling
#########################################################

# Set up graph for Dynamics and param histograms

fig4,ax4 = plt.subplots(2,3,figsize=[10,7])
fig4.suptitle('Abiotic HOOH Model Output')
#graph dynamics of data and model (best model line), and posterior guesses (grey lines) of 0 and 400 respectively 
ax4[0,0].plot(df0.time,df0.abundance, marker='o',label = 'abiotic - 0 H ')
ax4[0,0].plot(mod0.time,mod0['H'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax4[0,0],a1,posteriors1,100)
ax4[1,0].plot(df4.time,df4.abundance, marker='o',label = 'abiotic - 400 H ')
ax4[1,0].plot(mod4.time,mod4['H'],c='r',lw=1.5,label=' model best fit')
plot_uncertainty(ax4[1,0],a4,posteriors4,100)

fig4.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30)

l1 = ax4[0,0].legend(loc = 'upper left')
l2 = ax4[1,0].legend(loc = 'upper left')
l1.draw_frame(False)
l2.draw_frame(False)


# plot histograms of params next to dynamics graphs
ax4[0,1].hist(posteriors1.Sh)
ax4[0,2].hist(posteriors1.deltah)
ax4[1,1].hist(posteriors4.Sh)
ax4[1,2].hist(posteriors4.deltah)

#set titles of sebplots
ax4[0,0].set_title('Model-Data Dynamics')
ax4[0,1].set_title('Sh')
ax4[0,2].set_title('deltah')

plt.show()


#fig4.savefig('../figures/abiotic_odelib')

##############
#graph parameters against one another 
###############

#graph set up

fig5,ax5 = plt.subplots(2,1,sharex=True, figsize=[8,5])
fig5.suptitle('deltah vs Sh ')
#graphing each assay's parameters against each other 
ax5[0].scatter(posteriors1.Sh,posteriors1.deltah)
ax5[1].scatter(posteriors4.Sh,posteriors4.deltah)

#ax5[0].set_xlabel('Frequency Sh')
ax5[1].set_xlabel('Frequency Sh')
ax5[0].set_ylabel('Frequency deltah')
ax5[1].set_ylabel('Frequency deltah')
#ax5[0].set_yscale('log')
plt.legend()
plt.show()

#fig5.savefig('../figures/abiotic_params')


########
#graphing logged parameter values
#######

fig6,ax6 = plt.subplots(2,2,sharex=True,figsize=[8,5])
fig6.suptitle('Log deltah vs Sh ')
ax6[0,0].scatter(posteriors1.iteration,np.log(posteriors1.Sh))
ax6[0,1].scatter(posteriors4.iteration,np.log(posteriors4.Sh))
ax6[1,0].scatter(posteriors1.iteration,np.log(posteriors1.deltah))
ax6[1,1].scatter(posteriors4.iteration,np.log(posteriors4.deltah))

ax6[1,1].set_xlabel('iteration')
ax6[0,0].set_title('0 HOOH')
ax6[0,1].set_title('400 HOOH ')


ax6[0,0].set_ylabel('log Sh')
ax6[1,0].set_ylabel('log deltah')
#ax5[0].set_yscale('log')

plt.show()




fig7,ax7 = plt.subplots(2,2,sharex=True,figsize=[8,5])
fig7.suptitle('logged param exploration ')
ax7[0,0].scatter(posteriors1.rsquared,np.log(posteriors1.Sh))
ax7[0,1].scatter(posteriors4.rsquared,np.log(posteriors4.Sh))
ax7[1,0].scatter(posteriors1.rsquared,np.log(posteriors1.deltah))
ax7[1,1].scatter(posteriors4.rsquared,np.log(posteriors4.deltah))

ax7[1,1].set_xlabel('rsquared')
ax7[0,0].set_title('0 HOOH')
ax7[0,1].set_title('400 HOOH ')


ax7[0,0].set_ylabel('log Sh')
ax7[1,0].set_ylabel('log deltah')
#ax5[0].set_yscale('log')
l7 = ax7[0,0].legend()
l7.draw_frame(False)
plt.show()





print('we done did it')




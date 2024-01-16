'''

name:   model_cat_muts_noC_vibrio.py

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

df_all = df_MHMnoC


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

df = df[df['id'] =='CatMutsNoC_growth'] #set id of assay we are fitting

df0 = df.loc[(df['treatment(nM)']==0)]  #assay 0 H 
df15 = df.loc[(df['treatment(nM)']==1500)] #assay of spiked HOOH

dfw = df0   #setting working df


Ss = df['strain'].unique()
nSs = Ss.shape[0]
colors = ['k','lawngreen','forestgreen','seagreen','dodgerblue','deepskyblue','b','violet','fushia','mediumpurple','yellow','orange']


for (s,ns) in zip( Ss,range(nSs)):
    
    ##############
    #setting working directory
    ##############

    #slicing data into abiotic, biotic, and Pro only dataframes

    
    df = dfw[dfw['strain']==s]  #slicing working df by strain (different mutants) 

    ## Reading in inits files for 0 and 400 models respectively
    inits0 = pd.read_csv("../data/inits/MHMnoC_inits0.csv")


    #####################################################
    #functions  for modeling and graphing model uncertainty 
    #####################################################

    pw = 1

    #setting param prior guesses and inititaing as an odelib param class in odelib
    k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.0002})
    k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.5})
    #setting state variiable  prior guess
    D0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+4})
    N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':1e+6})
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
        return [dDdt,dNdt]
    
    
    #####################################################
    
    #find closest time 
    def get_residuals(self):
        mod = self.integrate(predict_obs=True)
        res = (mod.abundance - self.df.abundance)   #this is not same species 
        mod['res'] = res
        return(mod)
    
    
    
    # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
    nits = 10000
    
    
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
    c0 = colors[ns]
    
    
    # Set up graph for Dynamics and param histograms
    
    fig1, (ax0,ax1)= plt.subplots(1,2,figsize = (10,6)) #fig creationg of 1 by 2
    fig1.suptitle(str(s) +' in 0 H Model') #setting main title of fig
    
    ####### fig config and naming 
    
    fig1.subplots_adjust(right=0.9, wspace = 0.45, hspace = 0.20)
    
    ax0.semilogy()
    ax0.set_title(str(s)+' dynamics ')
    ax1.set_title('Model residuals')
    
    ax0.set_xlabel('Time (hrs)')
    ax0.set_ylabel('Cells (ml-1)')
    ax1.set_ylabel('Data D value')
    ax1.set_xlabel('Residual')


    #graphing data from df to see 2 different biological reps represented
    
    ax0.errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['rep1'],color = 'lightcoral', marker='^', label = 'rep1')
    ax0.errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['rep2'],color = 'firebrick' ,marker='.', label = 'rep2')
    ax0.plot(df[df['organism']== 'D']['time'], df[df['organism']== 'D']['rep3'], color = 'lightsalmon', marker='o',label = 'rep3')
    ax0.errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['abundance'],yerr=df[df['organism']== 'D']['sigma'],color = c0, marker='d', label = 'MEAN '+str(s))
    
    ax0.plot(mod0.time,mod0['D'],c='r',lw=1.5,label=' model best fit')
    
    a0.plot_uncertainty(ax0,posteriors0,'D',1000)
    
    ax1.scatter(a0res['res'], a0res['abundance'],color = c0,label = '0H case')
    
    l1 = ax0.legend(loc = 'lower right')
    l1.draw_frame(False)
    
    
    
    fig1.savefig('../figures/MHMnoC_'+s+'0_dynamics')
    
    
    
    ###################
    fig2,ax2 = plt.subplots(1,3,sharex=True,figsize=[8,5]) 
    
    ax2[0].set_title('Model Dynamic output', fontsize = 16)
    ax2[1].set_title('D0', fontsize = 12)
    ax2[2].set_title('k2', fontsize = 12)
    ax2[0].semilogy()
    ax2[1].set_xlabel('Parameter Value', fontsize = 12)
    ax2[1].set_ylabel('Frequency', fontsize = 12)
    ax2[1].xaxis.set_label_coords(0.88, -0.2)
    ax2[0].set_xlabel('Time (days)', fontsize = 14)
    #make legends

    #shift fig subplots
    fig2.subplots_adjust(right=0.95, wspace = 0.45, left = 0.05, hspace = 0.20, bottom = 0.2)
    
    
    #graph data, model, and uncertainty 
    #ax2[0].plot(df[df['organism']== 'D']['time'], df[df['organism']== 'D']['abundance'], marker='>',label = str(s) +' MEAN in 0 H ')
    ax2[0].errorbar(df[df['organism']== 'D']['time'],df[df['organism']== 'D']['abundance'],yerr=df[df['organism']== 'D']['sigma'], color = c0,marker='^', label = 'MEAN '+str(s))
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
    
    fig2.savefig('../figures/'+s+'0H_odelib')
    
    
    
    #################################
    #graphing trace of Param values
    ##################################
    #crating and config of fig 3
    fig3,ax3 = plt.subplots(1,2,sharex=True,figsize=[8,5]) #make plot
    fig3.suptitle('Trace plots for Params in 0H') #set main title 
    fig3.subplots_adjust(right=0.90, wspace = 0.55, top = 0.90) #shift white space for better fig view
    fig3.supxlabel('Model Iteration') #set overall x title 
    
    ax3[0].set_ylabel('D0')
    ax3[1].set_ylabel('k2')
    
    #graphing iteration number vs parameter numbert logged 
    ax3[0].scatter(posteriors0.iteration,(posteriors0.D0),color = c0)
    ax3[1].scatter(posteriors0.iteration,(posteriors0.k2),color = c0)
    
    #print out plot
    fig3.savefig('../figures/MHMnoC_'+s+'0_TRACE')
    
    
    plt.show()
    

    pframe0 = pd.DataFrame(a0.get_parameters(),columns=a0.get_pnames())
    pframe0.to_csv("../data/inits/MHMnoC_inits0.csv")

'''
##################################
#abtioic spike df1500
##################################
df = df15
inits15 = pd.read_csv("../data/inits/MHMnoC_inits1500.csv")



#setting param prior guesses and inititaing as an odelib param class in odelib
deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':1500})

priors = {'deltah' :deltah_prior,'Sh' : Sh_prior,'H0' : H0_prior} #list of all priors to feed to odelib create

#setting H mean for odelib search 
H0_mean = inits15['H0'][0]


# nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS
nits = 10000


#####################################
# Create and Run model on 0 and 400 df
#####################################

a15 = get_model(df15) 

posteriors15 = a15.MCMC(chain_inits=inits15,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )

# run model with optimal params
mod15 = a15.integrate()

a15res = get_residuals(a15)  #is this using the best fit or just a first run???


#########################################################
# graphing df and models together
#########################################################
c15 = 'darkviolet'


# Set up graph for Dynamics and param histograms

fig4,ax4 = plt.subplots(1,3,figsize=[9,4]) #plot creation and config 
#set titles of subplots
fig4.suptitle('Abiotic 1500 HOOH Model Output') #full title config
fig4.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30) #shift white space for better fig view
ax4[0].set_title('HOOH Dynamics')
ax4[0].set_ylabel('HOOH Concentration nM/mL')
ax4[0].set_xlabel('Time (days)')
ax4[1].set_title('Sh')
ax4[1].set_ylabel('Frequency')
ax4[1].set_xlabel('Parameter Value')
ax4[2].set_title('deltah')
ax4[2].set_ylabel('Frequency')
ax4[2].set_xlabel('Parameter Value (Logged)')


#plot dynamics of data and model for 0 assay 
ax4[0].plot(df15.time,df15.abundance, marker='o',color = c15, label = 'H data ') #data of 0 H assay
ax4[0].plot(mod15.time,mod15['H'],c='k',lw=1.5,label=' model best fit') #best model fit of 0 H assay
a0.plot_uncertainty(ax4[0],posteriors15,'H',100)

# plot histograms of params next to dynamics graphs
ax4[1].hist(((posteriors15.Sh)), facecolor=c15) #graphing Sh of 0 H assay 
ax4[2].hist(((posteriors15.deltah)), facecolor=c15) #graphing deltah of 0 H assay 

#config legends
l4 = ax4[0].legend(loc = 'lower right')
l4.draw_frame(False)


fig4.savefig('../figures/MHMnoC_vibrio1500_dynamics')


#################################
#graphing logged parameter values
##################################
#crating and config of fig 3
fig5,ax5 = plt.subplots(1,2,sharex=True,figsize=[8,5]) #make plot
fig5.suptitle('Trace plots for Logged Params in 1500 nM H') #set main title 
fig5.subplots_adjust(right=0.90, wspace = 0.55, top = 0.90) #shift white space for better fig view
fig5.supxlabel('Model Iteration') #set overall x title 

ax5[0].set_ylabel('Log Sh')
ax5[1].set_ylabel('Log deltah')

#graphing iteration number vs parameter numbert logged 
ax5[0].scatter(posteriors15.iteration,np.log(posteriors15.Sh),color = c15)
ax5[1].scatter(posteriors15.iteration,np.log(posteriors15.deltah),color = c15)

#print out plot
fig5.savefig('../figures/MHMnoC_vibrio1500_TRACE')


###############

#making and confing of residuals plot
fig6,ax6 = plt.subplots(figsize=[8,5])
fig6.suptitle('Residuals vs Model Fit Value  - abiotic 1500nM')
fig6.supylabel('Model Value (H)')
fig6.supxlabel('Residual')
#config legends for data differentialtion 
l6 = ax6.legend()
l6.draw_frame(False)

#plotting residual function output residual and abundance columns 
ax6.scatter(a15res['res'], a15res['abundance'],label = '1500nM H', color = c0) #where )

#how to get residuals from all posterior runs not just best???

#print out plot
fig6.savefig('../figures/MHMnoC_abiotic_1500_residuals')


###################



##############################
#vibrio modeling
##############################


#HOOH dynamics 
fig5,ax5 = plt.subplots(1,3,figsize=[10,5])
fig5.suptitle('HOOH parmaters ', fontsize = 14)
ax5[0].set_title('HOOH Model Dynamic output', fontsize = 14)
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


ax5[0].plot(df15[df15['organism']=='H']['time'], df15[df15['organism']=='H']['abundance'], marker='o',label = 'H data ')
ax5[0].plot(mod15.time,mod15['H'],color='r',lw=1.5,label=' Model H best fit')
a15.plot_uncertainty(ax5[0],posteriors15,'H',100)


# plot histograms of parameter search results 
ax5[1].hist(posteriors15.H0)
ax5[2].hist(posteriors15.phi)


#show full graph 

fig5.savefig('../figures/vibrio_odelib15_Hparams')


plt.show()


pframe15 = pd.DataFrame(a15.get_parameters(),columns=a15.get_pnames())
pframe15.to_csv("../data/inits/MHMnoC_inits1500.csv")
'''

plt.show()

# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print('\n Done my guy \n')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')



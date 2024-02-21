'''
viz_Pro_help_and_light.py

visualise different buffer proiductions of H from Morris 2013 (miliQ_bufferH_2013MorrisSIfig1b) 
  HOOH production in 10 mM Buffers plus seawater media incubated  in lightexposed milli-Q water (i.e., without seawater solutes). rep1-3 are bio1 rep4-6 are bio2

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   
import scipy 
import ODElib
import random as rd
import sys


#######functions

def get_model(df):
    M = ODElib.ModelFramework(ODE=mono_4H,
                              parameter_names=['k1','k2','kdam','phi','Sh','P0','N0','H0'],
                              state_names = snames,
                              dataframe=df,
                              k1 = k1_prior.copy(),
                              k2 = k2_prior.copy(),
                              kdam = kdam_prior.copy(),
                              phi = phi_prior.copy(),
                              Sh = Sh_prior.copy(),
                              P0 = P0_prior.copy(),
                              N0  = N0_prior.copy(),
                              H0  = H0_prior.copy(),
                              t_steps=1000,
                              P = P0_mean,
                            N = N0_mean,
                                    H = H0_mean,
                                        )
    return M

def mono_4H(y,t,params): #no kdam or phi here (or make 0)
    k1,k2, kdam, phi, Sh = params[0], params[1], params[2],params[3], params[4]
    P,N,H = max(y[0],0),max(y[1],0),y[2]
    ksp=k2/k1 #calculating model param ks in loop but k1 and k2 are fed separately by odelib
    dPdt = (k2 * N /( (ksp) + N) )*P - kdam*P*H    
    dNdt =  - (k2 * N /( (ksp) + N) )*P
    dHdt = Sh- phi*P*H
    return [dPdt,dNdt,dHdt]


def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)

###############################

#  Importing the data frame

###############################

 #UH18301 Pro in 3.75 mM HEPES or Taps buffer. High light (24uC in a Sunbox -  noon maximum of about 250 quanta m) 



df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = 'Pro_help_and_light_Morris2011_5', header  = 0)


df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_all = df_all.fillna(0)



df_all['rep1'] = df_all['rep1'].astype('float')
df_all['rep2']  = df_all['rep2'].astype('float')
df_all['rep3']  = df_all['rep3'].astype('float')
df = df_all 

#logging data for latter graphing 
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])

#####

#raw avgs reps

df['abundance'] =  np.nanmean(np.r_[[df[i] for i in ['rep1','rep2','rep3']]],axis=0)
df['sigma'] = np.nanstd(np.r_[[df[i] for i in ['rep1','rep2','rep3']]],axis=0)


#log avgs and stdvs

df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['log1','log2','log3']]],axis=0)
df['log_sigma'] =  np.nanstd(np.r_[[df[i] for i in ['log1','log2','log3']]],axis=0)

df['log_sigma'] = 0.2
df.loc[df['organism'] == 'H', 'log_sigma'] = 0.08

##############################

#   data arrays

###############################


orgs = df['Strain'].unique()
norgs = orgs.shape[0]


treats = df['Treatment'].unique()
ntreats = treats.shape[0]

lights = df['Light'].unique()
nlightss = lights.shape[0]
##############################

#    graphing the data 

##############################

for s,ns in zip(orgs,range(norgs)):
    for t,nt in zip(treats,range(ntreats)):
        dfw = df_all[(df_all['Strain'] == s) & (df_all['Treatment']==t)]
        dfw.loc[dfw['organism'] == 'S', 'organism'] = 'P'   #so Syn straing can be run here
        df_low =  dfw[(dfw['Light']== 'Low')]
        df_high = dfw[(dfw['Light']== 'High')]
        
    #df_ax = dfw[~(dfw['Treatment']== 'EZ55')]
    #df_help = dfw[(dfw['Treatment']== 'EZ55')]
        #print(dfw)

        fig1,(ax1)= plt.subplots(2,2, figsize = (11,8))
        fig1.suptitle(str(s) + ' ' + str(t), size = 22)
        fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
        ax1[0,0].set_ylabel('Pro (cells/ml)')
        ax1[0,0].set_xlabel('Time' )
        ax1[0,0].semilogy()
        ax1[1,0].set_ylabel('Pro (cells/ml)')
        ax1[1,0].set_xlabel('Time' )
        ax1[1,0].semilogy()
        ax1[1,1].set_ylabel('STDV')
        ax1[1,1].set_xlabel('Mean' )
        ax1[0,1].set_ylabel('STDV')
        ax1[0,1].set_xlabel('Mean' )

        fig2,(ax2) = plt.subplots(2,2,figsize = (11,8))
        fig2.suptitle(' Log ' +  str(s) + ' ' + str(t),size = 22)
        fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
        ax2[0,0].set_ylabel('Pro (cells/ml)')
        ax2[0,0].set_xlabel('Time' )
        ax2[0,0].semilogy()
        ax2[1,0].set_ylabel('Pro (cells/ml)')
        ax2[1,0].set_xlabel('Time' )
        ax2[1,0].semilogy()
        ax2[1,1].set_ylabel('Log STDV')
        ax2[1,1].set_xlabel('Log Mean' )
        ax2[0,1].set_ylabel('Log STDV')
        ax2[0,1].set_xlabel('Log Mean' )

        ax1[0,0].errorbar(df_low.time,df_low.abundance, yerr=df_low.sigma, marker = 'd', c='brown',label =  'raw mean')
        ax1[0,0].scatter(df_low.time,df_low.rep1, c='c', label = 'rep1')
        ax1[0,0].scatter(df_low.time,df_low.rep2, c='b', label = 'rep2')
        ax1[0,0].scatter(df_low.time,df_low.rep3, c='purple',label = ' rep3')
        ax1[0,0].text(1.2,0.5,'Low Light',horizontalalignment='center', verticalalignment='center', transform=ax1[0,1].transAxes)
        l1 = ax1[0,0].legend(loc = 'lower right')
        l1.draw_frame(False)
        ax1[0,1].scatter(df_low.abundance, df_low.sigma, c = 'purple')
    
        ax1[1,0].errorbar(df_high.time,df_high.abundance, yerr=df_high.sigma, marker = 'd', c='brown',label =  'raw mean')
        ax1[1,0].scatter(df_high.time,df_high.rep1, c='c', label = 'rep1')
        ax1[1,0].scatter(df_high.time,df_high.rep2, c='b', label = 'rep2')
        ax1[1,0].scatter(df_high.time,df_high.rep3, c='purple',label = ' rep3')
        ax1[1,1].text(1.2,0.5,'High Light',horizontalalignment='center', verticalalignment='center', transform=ax1[1,1].transAxes)
        #l1 = ax1[1,0].legend(loc = 'upper left')
        #l1.draw_frame(False)
        ax1[1,1].scatter(df_high.abundance,df_high.sigma, c = 'purple')
        
        #logged graph
        ax2[0,0].errorbar(df_low.time,df_low.log_abundance, yerr=df_low.log_sigma, marker = 'd', c='brown',label =  'log mean')
        ax2[0,0].scatter(df_low.time,df_low.log1, c='c', label = 'log1')
        ax2[0,0].scatter(df_low.time,df_low.log2, c='b', label = 'log2')
        ax2[0,0].scatter(df_low.time,df_low.log3, c='purple',label = ' log3')
        ax2[0,0].text(1.2,0.5,'Low Light',horizontalalignment='center', verticalalignment='center', transform=ax2[0,1].transAxes)
        l2 = ax2[0,0].legend(loc = 'lower right')
        l2.draw_frame(False)
        ax2[0,1].scatter(df_low.log_abundance, df_low.log_sigma, c = 'purple')
    
        ax2[1,0].errorbar(df_high.time,df_high.log_abundance, yerr=df_high.log_sigma, marker = 'd', c='brown',label =  'log mean')
        ax2[1,0].scatter(df_high.time,df_high.log1, c='c', label = 'log1')
        ax2[1,0].scatter(df_high.time,df_high.log2, c='b', label = 'log2')
        ax2[1,0].scatter(df_high.time,df_high.log3, c='purple',label = ' log3')
        ax2[1,1].text(1.2,0.5,'High Light',horizontalalignment='center', verticalalignment='center', transform=ax2[1,1].transAxes)
        #l3 = ax2[1,0].legend(loc = 'upper left')
        #l3.draw_frame(False)
        ax2[1,1].scatter(df_high.log_abundance,df_high.log_sigma, c = 'purple')
        
        fig1.savefig('../figures/'+str(s)+'raw'+str(t)+'.png')
        fig2.savefig('../figures/'+str(s)+'log'+str(t)+'.png')
        
        
        
    #model Axenic dynamics
    
        #if t == 'EZ55' :
            
        df4 = df_high
        if df4.shape[0] == 0:
            print('empty df in'+str(t)+str(s))
            continue
        if df4.shape[0] > 0:
            pass
#reading in csv file with inititla guesses for all parameter values ( SH, deltah, H0)
        inits4 = pd.read_csv('../data/inits/hepes_'+s+'_inits4.csv')

#setting how many MCMC chains you will run 
        nits = 1000 # nits - INCREASE FOR MORE BELL CURVEY LOOKING HISTS of params

# state variable names
        snames = ['P','N','H'] #order must match all further model mentions (same fro params) 

# define priors for parameters
        pw = 1   #sigma for param search


#setting param prior guesses and inititaing as an odelib param class in odelib
        k1_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':20000})
        k2_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.3})
        kdam_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.2})
        phi_prior = ODElib.parameter(stats_gen=scipy.stats.lognorm,hyperparameters={'s':pw,'scale':0.02})
        Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw,'scale':2})
#setting state variiable  prior guess
        P0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+4})
        N0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':2e+7})
        H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, hyperparameters={'s':pw/1,'scale':200})
#pw/10 for state variable initial conditions (P0, H0, N0) bc we theoretically have a better handle on thier values. (not completely holding constant like Qnp but not as loose as params either)

#still not sure what part of fitting algor this is used for
        P0_mean = inits4['P0'][0]
        N0_mean = inits4['N0'][0]
        #H0_mean = inits4['H0'][0]
        H0_mean = 200 #from asssay 
        #####################################################
        #functions  for modeling and graphing model uncertainty 
        #####################################################


#df0.loc[:,'log_abundance'] = np.log(10**df0.log_abundance)

# get_models
        a4 = get_model(df4) 

#broken here!!!!!!!!!!
# do fitting
        posteriors4 = a4.MCMC(chain_inits=inits4,iterations_per_chain=nits,cpu_cores=1,print_report=True) #, )
#posteriors1 = a1.MetropolisHastings(chain_inits=inits0,iterations_per_chain=nits,burnin = 500,cpu_cores=1,static_parameters=set(['Qnp']))

# run model with optimal params
        mod4 = a4.integrate()

#####################################################
# graphing model vs data in 0 H and associated error
#####################################################

###### fig set up
        fig3, ax3 = plt.subplots(1,2,figsize = (9,5)) #fig creationg of 1 by 2
        fig3.suptitle(str(s)+' in '+str(t)+' Model') #setting main title of fig

        ####### fig config and naming 
        
        fig3.subplots_adjust(right=0.85, wspace = 0.50, hspace = 0.30)

        ax3[0].semilogy()
        ax3[1].semilogy()
        ax3[0].set_title(str(s)+' dynamics ', fontsize = 14)
        ax3[1].set_title('HOOH dynamics', fontsize = 14)
        
        ax3[0].set_xlabel('Time (days)')
        ax3[0].set_ylabel('Cell concentration')
        ax3[1].set_ylabel('HOOH concentration')
        ax3[1].set_xlabel('Time (days)')
        
        

        #graphing data from df to see 2 different biological reps represented

        ax3[0].errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['abundance'],yerr=df4[df4['organism']=='P']['sigma'],c = 'g', marker='o', label = 'Mean P')
        ax3[0].plot(mod4.time,mod4['P'],color ='r',lw=1.5,label=' P model best fit')
        a4.plot_uncertainty(ax3[0],posteriors4,'P',100)

        ax3[1].errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['sigma'],c = 'goldenrod', marker='o', label = 'Mean H')
        ax3[1].plot(mod4.time,mod4['H'],color ='r',lw=2.0,label=' H model best fit')
        a4.plot_uncertainty(ax3[1],posteriors4,'H',100)
        l3 = ax3[0].legend(loc = 'lower right')
        l3.draw_frame(False)

        #ax1.scatter(a0res['res'], a0res['abundance'],label = '0H case')
        #printing off graph
        plt.show()

        fig3.savefig('../figures/pro_data_400')


        #########################################################
        #graphing P model vs data and params histograms 
        #########################################################

        # set up graph
        fig4,ax4 = plt.subplots(1,3,figsize=[10,4])
        #set titles and config graph 
        fig4.suptitle(str(s)+' Monoculture parameters in '+str(t))
        ax4[0].set_title(str(s)+' dyanmics')
        ax4[1].set_title('P0')
        ax4[2].set_title('kdam')

        ax4[0].semilogy()
        ax4[0].set_ylabel('Cell concentration')
        ax4[0].set_xlabel('Time (Days)')


        ax4[1].set_xlabel('Parameter Value', fontsize = 12)
        ax4[1].set_ylabel('Frequency', fontsize = 12)
        #make legends

        #shift fig subplots
        fig4.subplots_adjust(right=0.90, wspace = 0.45, hspace = 0.30)


        #graph data, model, and uncertainty 
        ax4[0].errorbar(df4[df4['organism']=='P']['time'],df4[df4['organism']=='P']['abundance'],yerr=df4[df4['organism']=='P']['sigma'],c = 'g', marker='o', label = 'Mean P')
        ax4[0].plot(mod4.time,mod4['P'],color='r',lw=1.5,label=' Model P best fit')
        a4.plot_uncertainty(ax4[0],posteriors4,'P',100)

        # plot histograms of parameter search results 
        ax4[1].hist(posteriors4.P0, color =  'g')
        ax4[2].hist(posteriors4.kdam,color =  'g')
        
        l4 = ax4[0].legend(loc = 'upper left')
        l4.draw_frame(False)
        #show full graph 
        plt.show()
        fig4.savefig('../figures/'+str(s)+'_hepes_model_Pparams')

        '''
        #########################################################
        #graphing H model vs data and params histograms 
        #########################################################

        #HOOH dynamics 
        fig5,ax5 = plt.subplots(1,3,figsize=[10,4])
        fig5.suptitle('HOOH parmaters in HEPES')
        ax5[0].set_title('HOOH dynamics')
        ax5[1].set_title('H0')
        ax5[2].set_title('phi')

        ax5[0].set_ylabel('HOOH concentration')
        ax5[0].set_xlabel('Time (Days)')
        fig5.subplots_adjust(right=0.85, wspace = 0.45, hspace = 0.20)

        ax5[1].set_xlabel('Parameter Value', fontsize = 12)
        ax5[1].set_ylabel('Frequency', fontsize = 12)

        ax5[0].errorbar(df4[df4['organism']=='H']['time'],df4[df4['organism']=='H']['abundance'],yerr=df4[df4['organism']=='H']['sigma'],c = 'g', marker='o', label = 'Mean H')
        ax5[0].plot(mod4.time,mod4['H'],color='r',lw=1.5,label=' Model H best fit')
        a4.plot_uncertainty(ax5[0],posteriors4,'H',100)


        # plot histograms of parameter search results 
        ax5[1].hist(posteriors4.H0,color =  'goldenrod')
        ax5[2].hist(posteriors4.phi, color =  'goldenrod')


        l5 = ax5[0].legend(loc = 'upper left')
        l5.draw_frame(False)

        #show full graph 
        plt.show()
        fig5.savefig('../figures/'+str(s)+'_hepes_model_Hparams')
        '''
 
        pframe = pd.DataFrame(a4.get_parameters(),columns=a4.get_pnames())
        pframe.to_csv('../data/inits/hepes_'+str(s)+'_inits4.csv')

    else: 
        print(t)
    
    

    


    

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.show()
#fig1.savefig('../figures/Pro_help_with_light_raw.png')
#fig2.savefig('../figures/logDynamics_Pro_help_with_light.png')




# 'program finished' flag

print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')
print('done with Pro light and help assays')
print('\n ~~~****~~~****~~~ \n')

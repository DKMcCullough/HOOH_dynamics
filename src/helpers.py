from scipy.stats.distributions import chi2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy

plt.rcParams["font.family"] = "DeJavu Serif"
pd.Series.iteritems = pd.Series.items

def get_data(identifier,sheet='BCC_1-31-dataset'):
    df_all = pd.read_excel("../data/ROS_data_MEGA.xlsx",sheet_name = sheet, header = 1)
    df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    df_all = df_all.rename({'time(day)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
    df_a = df_all.loc[df_all['assay'].str.contains(identifier, case=False)].copy()
    if identifier == 'coculture':
        df_a = df_all.loc[~df_all['assay'].str.contains(identifier, case=False)].copy()
    df = summary_stats(df_a)
    df['log_sigma'] = 0.056225 # default log sigma
    return df

def get_uncertainty(df,alpha=0.01):
    ddf = df.shape[0] - 1
    return np.mean(np.std(df.log3 - df.lavg1)*np.sqrt((ddf)/chi2.ppf(alpha/2,df=ddf)))

def summary_stats(df):
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
    return df

def plot_uncertainty(df,sigma0=0.2):
    figa,axa = plt.subplots(1,2,figsize=[9,4])
    axa[0].plot(df.time,df.log1,marker='o',label='rep 1')
    axa[0].plot(df.time,df.log3,marker='o',label='rep 2')
    axa[1].plot(df.lavg1,df.stdlog1,marker='o')
    axa[1].axhline(y=sigma0,c='r',label='99% CI upper limit')
    axa[1].axhline(y=np.mean(df.stdlog1),c='r',ls='--',label='mean across samples')
    axa[0].set_xlabel('Time (days)')
    axa[1].set_xlabel('Mean of logged replicates')
    axa[0].set_ylabel('Log abundance')
    axa[1].set_ylabel('Standard deviation of logged replicates')
    l1 = axa[0].legend(loc='lower right')
    l2 = axa[1].legend()
    l1.draw_frame(False)
    l2.draw_frame(False)
    for (l,ax) in zip('ab',axa):
        ax.text(0.07,0.9,l,ha='center',va='center',color='k',transform=ax.transAxes)
    return figa,axa



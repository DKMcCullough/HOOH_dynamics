#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:13:24 2023
Created on Tue Nov 15 23:50:19 2022

DDname:phi_vs_kdam.py

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/HOOH_dynamics/src

author: DKM


goal: icompare phi and kdam values 
"""



#SynWH7803 (Vol 28 and 52) is KatG possitive Syn + 
#SynWH8102 (Vol 53) is KatG negative  Syn - 
#Pro = MIT9215

import helpers as hp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

####################################################
# optimal params
#####################################################

# spike
df_Pro_1 =  pd.read_csv('../data/inits/pro_MIT9215_inits4_1.csv') #Pro MIT9215
df_Pro_2 = pd.read_csv('../data/inits/pro_MIT9215_inits4_2.csv') #Pro MIT9215
df_Pros = pd.concat([df_Pro_1, df_Pro_2]).groupby(level=0).mean()   #avg of Pro MIT9215
df_Syn_28 = pd.read_csv('../data/inits/syn_vol28_inits4.csv') #Syn WH7803
df_Syn_52 = pd.read_csv('../data/inits/syn_vol52_inits4.csv') #Syn WH7803
df_Syn_WH7803 = pd.concat([df_Syn_28 , df_Syn_52]).groupby(level=0).mean()   #avg of Syn WH7803
df_Syn_53 = pd.read_csv('../data/inits/syn_vol53_inits4.csv') #Syn CC9605
df_Het_57 = pd.read_csv('../data/inits/Het_57_inits4.csv') #Micromonas commoda
df_Het_58 = pd.read_csv('../data/inits/Het_58_inits4.csv') #Micromonas pusilla
df_Het_59 = pd.read_csv('../data/inits/Het_59_inits4.csv') #Ostreococcus lucimarinus
df_Het_60 = pd.read_csv('../data/inits/Het_60_inits4.csv') #Ostreococcus tauri
df_Het_55 = pd.read_csv('../data/inits/Het_55_inits_300-spike.csv') #Alteromonas macleodii  in 300nM 

df_spike = pd.DataFrame({'names' : ('$Synechococcus$ WH7803' ,'$Synechococcus$ CC9605', '$Prochlorococcus$ MIT9215' ,'$Micromonas$ $commoda$','$Micromonas$ $pusilla$','$Ostreococcus$ $lucimarinus$','$Ostreococcus$ $tauri$', '$Alteromonas$ $macleodii$'),
                                'phis' : [df_Syn_WH7803.phi,df_Syn_53.phi,df_Pros.phi,df_Het_57.phi, df_Het_58.phi,df_Het_59.phi,df_Het_60.phi,df_Het_55.phi],
                                'kdams' : [df_Syn_WH7803.kdam,df_Syn_53.kdam,df_Pros.kdam,df_Het_57.kdam, df_Het_58.kdam,df_Het_59.kdam,df_Het_60.kdam ,df_Het_55.kdam], 
                                'markers' : ['o','o','D','v','v','v','v','d'],
                                'colors' : ['dodgerblue','steelblue','lightgreen','violet','blueviolet','mediumorchid','mediumpurple','k']}, 
                                columns=['names','phis', 'kdams','markers','colors'])

# no spike
df_Pros =  pd.read_csv('../data/inits/pro_MIT9215_inits0_1.csv') #Pro MIT9215
df_Syn_28 = pd.read_csv('../data/inits/syn_vol28_inits0.csv') #Syn WH7803
df_Syn_52 = pd.read_csv('../data/inits/syn_vol52_inits0.csv') #Syn WH7803
df_Syn_WH7803 = pd.concat([df_Syn_28 , df_Syn_52]).groupby(level=0).mean()   #avg of Syn WH7803
df_Syn_53 = pd.read_csv('../data/inits/syn_vol53_inits0.csv') #Syn CC9605
df_Het_57 = pd.read_csv('../data/inits/Het_57nospike_inits4.csv') #Micromonas commoda
df_Het_58 = pd.read_csv('../data/inits/Het_58nospike_inits4.csv') #Micromonas pusilla
df_Het_59 = pd.read_csv('../data/inits/Het_59nospike_inits4.csv') #Ostreococcus lucimarinus
df_Het_60 = pd.read_csv('../data/inits/Het_60nospike_inits4.csv') #Ostreococcus tauri
df_Het_55 = pd.read_csv('../data/inits/Het_55_inits_0-spike.csv') #Alteromonas macleodii  in 300nM 

df_nospike = pd.DataFrame({'names' : ('$Synechococcus$ WH7803' ,'$Synechococcus$ CC9605', '$Prochlorococcus$ MIT9215' ,'$Micromonas$ $commoda$','$Micromonas$ $pusilla$','$Ostreococcus$ $lucimarinus$','$Ostreococcus$ $tauri$', '$Alteromonas$ $macleodii$'),
                                'phis' : [df_Syn_WH7803.phi,df_Syn_53.phi,df_Pros.phi,df_Het_57.phi, df_Het_58.phi,df_Het_59.phi,df_Het_60.phi,df_Het_55.phi],
                                'kdams' : [df_Syn_WH7803.kdam,df_Syn_53.kdam,df_Pros.kdam,df_Het_57.kdam, df_Het_58.kdam,df_Het_59.kdam,df_Het_60.kdam ,df_Het_55.kdam], 
                                'markers' : ['o','o','D','v','v','v','v','d'],
                                'colors' : ['dodgerblue','steelblue','lightgreen','violet','blueviolet','mediumorchid','mediumpurple','k']}, 
                                columns=['names','phis', 'kdams','markers','colors'])

####################################################
# posteriors
#####################################################

#spike
pos_Pro_1 =  pd.read_csv('../data/posteriors/pos_pro_spike.csv') #Pro MIT9215
pos_Pro_2 = pd.read_csv('../data/posteriors/pos_pro_spike2.csv') #Pro MIT9215
pos_Pros = pd.concat([pos_Pro_1])   #avg of Pro MIT9215
pos_Syn_28 = pd.read_csv('../data/posteriors/pos_syn_28spike.csv') #Syn WH7803
pos_Syn_52 = pd.read_csv('../data/posteriors/pos_syn_52spike.csv') #Syn WH7803
pos_Syn_WH7803 = pd.concat([df_Syn_28 , df_Syn_52]) # concatenate
pos_Syn_53 = pd.read_csv('../data/posteriors/pos_syn_53spike.csv') #Syn CC9605
pos_Het_57 = pd.read_csv('../data/posteriors/pos_het_57spike.csv') #Micromonas commoda
pos_Het_58 = pd.read_csv('../data/posteriors/pos_het_58spike.csv') #Micromonas pusilla
pos_Het_59 = pd.read_csv('../data/posteriors/pos_het_59spike.csv') #Ostreococcus lucimarinus
pos_Het_60 = pd.read_csv('../data/posteriors/pos_het_60spike.csv') #Ostreococcus tauri
pos_Het_55 = pd.read_csv('../data/posteriors/pos_EZ55_300-spike.csv') #Alteromonas macleodii  in 300nM 

pos_Syn_WH7803 = pos_Syn_28

spike_df = {
    df_spike['names'][2]: pos_Pros,
    df_spike['names'][0]: pos_Syn_WH7803,
    df_spike['names'][1]: pos_Syn_53,
    df_spike['names'][3]: pos_Het_57,
    df_spike['names'][4]: pos_Het_58,
    df_spike['names'][5]: pos_Het_59,
    df_spike['names'][6]: pos_Het_60,
    df_spike['names'][7]: pos_Het_55,
}

#no spike
pos_Pros =  pd.read_csv('../data/posteriors/pos_pro_no_spike.csv') #Pro MIT9215
pos_Syn_28 = pd.read_csv('../data/posteriors/pos_syn_28nospike.csv') #Syn WH7803
pos_Syn_52 = pd.read_csv('../data/posteriors/pos_syn_52nospike.csv') #Syn WH7803
pos_Syn_WH7803 = pd.concat([df_Syn_28 , df_Syn_52]) # concatenate
pos_Syn_53 = pd.read_csv('../data/posteriors/pos_syn_53nospike.csv') #Syn CC9605
pos_Het_57 = pd.read_csv('../data/posteriors/pos_het_57nospike.csv') #Micromonas commoda
pos_Het_58 = pd.read_csv('../data/posteriors/pos_het_58nospike.csv') #Micromonas pusilla
pos_Het_59 = pd.read_csv('../data/posteriors/pos_het_59nospike.csv') #Ostreococcus lucimarinus
pos_Het_60 = pd.read_csv('../data/posteriors/pos_het_60nospike.csv') #Ostreococcus tauri
pos_Het_55 = pd.read_csv('../data/posteriors/pos_EZ55_0-spike.csv') #Alteromonas macleodii  in 300nM 

pos_Syn_WH7803 = pos_Syn_28

nospike_df = {
    df_spike['names'][2]: pos_Pros,
    df_spike['names'][0]: pos_Syn_WH7803,
    df_spike['names'][1]: pos_Syn_53,
    df_spike['names'][3]: pos_Het_57,
    df_spike['names'][4]: pos_Het_58,
    df_spike['names'][5]: pos_Het_59,
    df_spike['names'][6]: pos_Het_60,
    df_spike['names'][7]: pos_Het_55,
}

####################################################
#monoculture of all  Pro, Syn, and Het Kdam and Phis
#####################################################

fig3, (ax1,ax2)= plt.subplots(1,2,figsize = (15,6))
ax1.set_xlabel('Detoxification rate x10$^{-6}$ ($\phi_{det,i})$',fontsize = 16)
ax1.set_ylabel('Damage rate ($\kappa_{dam,i}$)',fontsize = 16)
ax1.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

ax2.set_xlabel('Detoxification rate x10$^{-6}$ ($\phi_{det,i})$',fontsize = 16)
ax2.set_ylabel('Damage rate ($\kappa_{dam,i}$)',fontsize = 16)
ax2.tick_params(axis = 'both', which = 'both', length = 5, labelsize = 16)

A_COL = "phi"
B_COL = "kdam"

nospike_rows, spike_rows = [],[]

for label, df in spike_df.items():
    a_vals = df[A_COL].to_numpy()
    b_vals = df[B_COL].to_numpy()

    a_q25, a_q50, a_q75 = np.percentile(a_vals, [0, 50, 100])*1e+6
    b_q25, b_q50, b_q75 = np.percentile(b_vals, [0, 50, 100])

    spike_rows.append({
        "label": label,
        "a_med": a_q50,
        "b_med": b_q50,
        "a_err_low":  a_q50 - a_q25,
        "a_err_high": a_q75 - a_q50,
        "b_err_low":  b_q50 - b_q25,
        "b_err_high": b_q75 - b_q50,
    })

for label, df in nospike_df.items():
    a_vals = df[A_COL].to_numpy()
    b_vals = df[B_COL].to_numpy()

    a_q25, a_q50, a_q75 = np.percentile(a_vals, [0, 50, 100])*1e+6
    b_q25, b_q50, b_q75 = np.percentile(b_vals, [0, 50, 100])

    nospike_rows.append({
        "label": label,
        "a_med": a_q50,
        "b_med": b_q50,
        "a_err_low":  a_q50 - a_q25,
        "a_err_high": a_q75 - a_q50,
        "b_err_low":  b_q50 - b_q25,
        "b_err_high": b_q75 - b_q50,
    })


spike_summary = pd.DataFrame(spike_rows)
nospike_summary = pd.DataFrame(nospike_rows)

# Scatter via seaborn
sns.scatterplot(
    data=nospike_summary,
    x="a_med", y="b_med",
    ax=ax1,
    s=40
)
sns.scatterplot(
    data=spike_summary,
    x="a_med", y="b_med",
    ax=ax2,
    s=40
)

# 2D error bars via matplotlib
ax1.errorbar(
    nospike_summary["a_med"],
    nospike_summary["b_med"],
    xerr=[nospike_summary["a_err_low"], nospike_summary["a_err_high"]],
    yerr=[nospike_summary["b_err_low"], nospike_summary["b_err_high"]],
    fmt="none",        # don't add extra markers
    ecolor="gray",
    elinewidth=1,
    capsize=3,
    alpha=0.8
)
ax2.errorbar(
    spike_summary["a_med"],
    spike_summary["b_med"],
    xerr=[spike_summary["a_err_low"], spike_summary["a_err_high"]],
    yerr=[spike_summary["b_err_low"], spike_summary["b_err_high"]],
    fmt="none",        # don't add extra markers
    ecolor="gray",
    elinewidth=1,
    capsize=3,
    alpha=0.8
)


for i, r in df_nospike.iterrows():
    ax1.plot(r['phis']*1e+6, r['kdams'], 'o',c = r['colors'], marker = r['markers'],markersize=11, linewidth=0.1, label=r['names'])
for i, r in df_spike.iterrows():
    ax2.plot(r['phis']*1e+6, r['kdams'], 'o',c = r['colors'], marker = r['markers'],markersize=11, linewidth=0.1, label=r['names'])


l3 = ax1.legend(loc = 'best', fontsize = 12)

l3.draw_frame(False)
fig3.subplots_adjust(bottom=0.15)

ax1.semilogx()
ax1.semilogy()

ax2.semilogx()
ax2.semilogy()

for a in (ax1,ax2):
    a.set_xlim([1e-5,1e+2])
    a.set_ylim([1e-5,1])

fig3.savefig('../figures/tradeoffs',bbox_inches='tight')
fig3.savefig('../figures/figure8.tiff',bbox_inches='tight',dpi=300,format='tiff')

f4,ax4 = plt.subplots(8,2,figsize=[15,15])

for ax in ax4.flat:
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# re-enable y ticks + labels on left column only
for ax in ax4.flat:
    ax.tick_params(
        axis="y",
        left=True,
        labelleft=True
    )

# re-enable x ticks + labels on bottom row only
for ax in ax4[-1, :]:
    ax.tick_params(
        axis="x",
        bottom=True,
        labelbottom=True
    )

for (n,s,a,b) in zip(nospike_df.keys(),spike_df.keys(),ax4[:,0],ax4[:,1]):
    sns.kdeplot(np.log10(nospike_df[n].kdam),label='Spike Free',ax=a,linewidth=3)
    sns.kdeplot(np.log10(spike_df[n].kdam),label='Spiked',ax=a,linewidth=3)
    sns.kdeplot(np.log10(nospike_df[n].phi),label='Spike Free',ax=b,linewidth=3)
    sns.kdeplot(np.log10(spike_df[n].phi),label='No Spike',ax=b,linewidth=3)
    a.set_xlim([-8,-2])
    b.set_xlim([-8,-4])

ax4[0,0].text(0.3, 0.8, list(spike_df.keys())[0] ,transform=ax4[0,0].transAxes,va="center",ha="center",fontsize=14)
ax4[1,0].text(0.75, 0.8, list(spike_df.keys())[1] ,transform=ax4[1,0].transAxes,va="center",ha="center",fontsize=14)
ax4[2,0].text(0.75, 0.8, list(spike_df.keys())[2] ,transform=ax4[2,0].transAxes,va="center",ha="center",fontsize=14)
ax4[3,0].text(0.75, 0.8, list(spike_df.keys())[3] ,transform=ax4[3,0].transAxes,va="center",ha="center",fontsize=14)
ax4[4,0].text(0.75, 0.8, list(spike_df.keys())[4] ,transform=ax4[4,0].transAxes,va="center",ha="center",fontsize=14)
ax4[5,0].text(0.75, 0.8, list(spike_df.keys())[5] ,transform=ax4[5,0].transAxes,va="center",ha="center",fontsize=14)
ax4[6,0].text(0.75, 0.8, list(spike_df.keys())[6] ,transform=ax4[6,0].transAxes,va="center",ha="center",fontsize=14)
ax4[7,0].text(0.75, 0.8, list(spike_df.keys())[7] ,transform=ax4[7,0].transAxes,va="center",ha="center",fontsize=14)

ax4[3,0].text(-0.1, 0.5, 'Density',rotation=90,transform=ax4[3,0].transAxes,va="center",ha="center",fontsize=14)
ax4[0,0].legend()
ax4[-1,0].set_xlabel(r'log10 Damage rage, $\kappa_{dam,i}$',fontsize=14)
ax4[-1,1].set_xlabel(r'log10 Detoxification rate, $\phi_{det,i}$',fontsize=14)

for a in ax4.flat:
    a.set_ylabel("")
    a.tick_params(axis='both', which='major', labelsize=14) # Change 'major' to 'minor' for minor ticks

f4.subplots_adjust(wspace = 0.05,hspace = 0.1)

f4.savefig('../figures/dists.png',bbox_inches='tight')

print('Done!!!')

'''
Hepes_production_of_hooh.py

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Using Morrisand Zinser 2013 Fig 1_a  data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Hepes_Morris_2013/scripts/
'''


from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




###############################

#  Importing the data frame

###############################

df_all = pd.read_csv("../data/Hepes_assay_2013.csv") #use relative paths 
df_all.rename(columns = {'time(days)':'time'},inplace = True) 
df_all.dropna()
#xls = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx')
#hepes_10_df = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx','10mM Hepes',delimiter = '.', header = 2, names = ( 'Time (days)', 'HOOH Concentration (\u03BCM)')
#10mM_HOOH_df = pd.read_csv('../data/data_fig1a_Morris_and_Zinser_2013', delimiter =',', header= 2, names =( 'Time (1/days)','HOOH concentration'))


HOOH_df = np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4', 'rep5', 'rep6']]],axis=0)


##############################

#   data arrays

###############################

HOOH_data = 10**np.array(HOOH_df)   #raised to 10 because ot the way data thief recorded the data and its scale

#print(HOOH_data)





##############################

#   time arrays

##############################


HOOH_times =np.array(df_all['time'])
#print(HOOH_times)





##############################

#    graphing the data 

##############################


plt.scatter(x = HOOH_times, y = HOOH_data, marker = 'x', s=50, c = 'r', label = 'Measured HOOH') 
plt.xlabel('Time (day $^{-1}$)', fontsize=16) 
plt.ylabel('HOOH concentration (\u03BCM)',fontsize=16)
plt.yscale('log')
plt.tick_params(labelsize=12)
#plt.show()

'''
Hs = HOOH_data
Ts = HOOH_times

Hsd = np.std(Hs)     #getting numpy to find the stdv of the data (but really will probs be 0.1 bc the data set doesn't have triplicates in it. 




####################################

#analytical solution

####################################

#initial values and creating time array

delta = 0.6
S_HOOH = 2.8
step = 0.05 #delta t
ndays = 7
times = np.linspace(0,ndays,int(ndays/step))

def f(t, S_HOOH, delta):
	H = (S_HOOH/delta)*(1-e**(-delta*t))
	return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


#plt.plot(times,Hs,c='g',marker='*',label='Analytical Solution')






#############################

#Euler's Integration

############################

HsEuler = np.array([]) 
H = 0.0
t0 = 0

for t in times: 
	HsEuler = np.append(HsEuler,H)
	dHdt = S_HOOH - delta*H
	H = H + dHdt*step
	
#plt.plot(times,HsEuler,c='k',label = "Euler's Aproximation")

plt.legend()

#plt.show()



####################################

#ODE int

####################################


from scipy.integrate import odeint

def HsODEint(H,t):
	dHdt = S_HOOH-delta*H
	#print(dHdt,t,H)
	return dHdt


ode_solutions = odeint(HsODEint,0,times)


#plt.plot(times,ode_solutions,c='purple', linestyle = ':', label = 'Odeint Approximation')
'''
plt.legend()

plt.show()




print('Done')


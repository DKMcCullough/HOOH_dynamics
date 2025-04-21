#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: functions_chap2.py
functions to be used for monoculture modeling of abioitc H202, biotic growth, and biotic h202 interactions
@author: dkm
"""



#Base parameter set 


Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 
Qns = 1#(20.0e-15*(1/(14.0))*1e+9) 

k1p =  0.02     #Pro alpha
k1s =  0.01    #Syn alpha 
k2p =  0.5    #Vmax     P
k2s =  0.5  #Vmax     S  #0.300748361
ksp = k2p/k1p
kss = k2s/k1s
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.005   #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phi = 1.7e-6    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002
Sh = 0
SN = 10000

params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

##############
#make param for heterotroph detox value and heterotroph to be unchanging with current N and also to be callable into leaky function
#############

#functions for model 

def leak(y,t,params):
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh = params[0], params[1], params[2], params[3],params[4],  params[5], params[6], params[7],params[8],params[9], params[10],params[11]
    P,S,N,H = y[0],y[1],y[2],y[3]
    dPdt = (k2p * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt =(k2s * N /( (kss) + N))*S - (ds *S) #- kdams*H*S      
    dNdt =  SN - ((k2p * N /( (ksp) + N) )*P* Qnp) - ((k2s * N /( (kss) + N))*S* Qns) - rho*N    
    dHdt = Sh - deltah*H  - phi*S*H  - (0.0004)*(1e5)*(H) #phi being S cell-specific detox rate
    return [dPdt,dSdt,dNdt,dHdt]

#extra detox 

def Pwins (params): 
    Nstar = ((ksp*dp )+(ksp*kdam))/((k2p*Qnp) - dp - kdam)
    Pstar = (SN - rho*Nstar)*((Nstar + ksp)/((k2p*Nstar)*Qnp))
    Sstar = 0
    Hstar = Sh/(deltah+phi*Pstar)  #do we need toassume H must be 0 for P to win?????
    return  Nstar, Pstar, Sstar, Hstar 



def Swins (params): 
    Nstar = (ds*kss)/((k2s*Qns)-ds)
    Sstar = (SN - rho*Nstar)*(((Nstar + kss)/(k2s*Nstar*Qns)))
    Pstar = 0
    Hstar = Sh/(deltah+phi*Pstar)  #do we need toassume H must be 0 for P to win?????
    return  Nstar, Pstar, Sstar, Hstar 



def Coexist (params): 
    Nstar = (kss*ds)/(k2s-ds)
    Pstar = ((SN-rho*Nstar)*(Nstar + ksp))/(k2p*Nstar*Qnp)
    Hstar = (((k2p*Nstar)/(Nstar + ksp))-(dp))*(1/kdam)
    Sstar = (Sh - deltah*Hstar)/(phi*Hstar)
    return  Nstar, Pstar, Sstar, Hstar 


def StarContour(Nstar, Pstar, Sstar, Hstar):
    #
    return
#N threshold 

#Nstarph = ((ksp*dp )+(ksp*kdam*Hstar))/((k2p*Qnp) - dp - (kdam*Hstar))

#h thresthold 
#vHline = ((deltah)/(Pstar*kdam)*((Nstarp+ksp)/(k2p*Nstarp*Pstar*Qnp)+(dp*Pstar)))


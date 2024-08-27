#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:21:06 2023

@author: jarrettvalenti
"""

#Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#economy Initial Values
T=5000#num years
K_0=135#67.79 #inital gwp
N_0=6838#initial population
A_0=3.8#initial tech factor
I_0=16.38#initial investment
Y_0=67.78
gamma=0.3
delta=0.1
Q_0=67.79
beta=0.0485

#climate initial vals pulled from DICE
psi_1=0.00267
psi_2=2
Tat_0=0.8
Tlo_0=0.0068
z_1=0.098
z_2=1.31
z_3=0.088
z_4=0.025
O_0= (1-(1/(1+psi_1*Tat_0 + psi_2*Tat_0**(2))))

sigma_0=0.549
sigma=.01#change to sigma
mu_0=.001*sigma_0
theta_1=sigma_0/2800
theta_2=2.8
V_0=theta_1*(mu_0**theta_2)

Mat_0=588 #atmosphere initial co2
Muo_0=1350 #upper ocean initial co2
Mlo_0=10010 #inital lower ocean co2
phi_at_at=0.912
phi_at_uo=0.088
phi_uo_at=0.0383
phi_uo_uo=0.9592
phi_uo_lo=0.0025
phi_lo_uo=0.0003375
phi_lo_lo=0.9996625

E_0=280 #historical value of ghg
E_1=830 #inital value of ghg
Eland_0=3.3#initial Eland co2
eps=.1 #emissions/Y

fex=0.7

#Shared resources initial vals
US_0=10596+31471+6959#reserves in exajoules as of 2020
us_var=131+151+174 #consumption in exajoules 2020
S_0= 77.5#around 17% of energy generated each year is renewable as of 2020
LA_0=4060000000#hectares of undeveloped land
NR_0=S_0+LA_0+US_0#sum
s_var=.05#growth of renewables
S_var= 40#current consumption of renewables in exajoules
l_var=0.0015#rate of aforestation
L_var=10000000#hectares consumed per year as of 2010-2020
dl=.13#ratio of econ attributed to land use

#Arrays to hold values
Y=[]
A=[]
K=[]
I=[]
N=[]
E=[]
Eland=[]
Q=[]
Mat=[]
Muo=[]
Mlo=[]
Tat = []
Tlo=[]
O=[]#damages
V=[]#abatement
F=[]#forcings
#shared resources arrays
NR=[]
S=[]
US=[]
LA=[]

for t in range(0,T):
    if(t==0):
        #nat resc
        nr=NR_0
        s=max(S_0*(1 + s_var) - S_var*(1 + dl),0)
        us=US_0
        la=LA_0
        #econ
        a= A_0
        k= K_0
        i= I_0
        n= N_0
        y= Y_0
        #clim
        o=O_0
        v=V_0
        sig = sigma_0
        mu= mu_0
        eland=Eland_0
        e=E_0
        mat=Mat_0
        muo=Muo_0
        mlo=Mlo_0
        f=5.33*np.log((mat)/Mat_0)/np.log(2) + fex
        tat=Tat_0
        tlo=Tlo_0
        q=Q_0
    else:
        #nat resc
        us = max(US_0 - us_var*(1+beta)**(t/1.5),0)
        if(us==0):#when out of resources demand will shift to renewables
            s = min(S_0*(1 + s_var)**(np.log(t)) - (S_var+us_var)*(1 + beta)**(25*np.log(t)),#75*s_var
                    1.2*(S_var+us_var)*(1 + beta)**(25*np.log(t)))
        else:
            s = S_0*(1 + s_var)**(25*np.log(t)) - S_var*(1 + beta)**(25*np.log(t))
        la = max((1 + l_var)*LA[t-1] - (L_var*(1+beta*dl)**(t)),0)
        nr= (s + la+ us)
        
        #econ
        a = np.log(3.8 + 0.079*t) +2.5
        i= 16.38*(1.01)**np.log((t))
        k= (nr/NR[t-1] - delta)*K[t-1] + (i)
        n= 6838 + 1.9*t
        y = (a*(k**gamma)*(np.log(n**(1-gamma))))-30
        
        #climate
        o=(1-(1/(1+psi_1*Tat[t-1]+psi_2*Tat[t-1]**(2))))
        sig = sigma_0+1.1**(sigma*t)
        mu = 0.001*sig
        eland = max(Eland[t-1] + 0.2*(1-la/LA[0]),0)
        e = sig*(1 - mu)*y*nr/NR[t-1] + eland - (s_var + l_var)*E[t-1] + 176
        mat = phi_at_at*Mat[t-1] + phi_uo_at*Muo[t-1]
        muo = phi_at_uo*Mat[t-1] + phi_uo_uo*Muo[t-1] + phi_lo_uo*Mlo[t-1]
        mlo = phi_uo_lo*Muo[t-1] + phi_lo_lo*Mlo[t-1]
        f = 5.33*(np.log(mat/Mat_0)/np.log(2)) + fex
        tat = Tat[t-1] + z_1*(f - z_2*Tat[t-1] - z_3*(Tat[t-1]-Tlo[t-1]))
        tlo = Tlo[t-1] + z_4*(Tat[t-1]-Tlo[t-1])
        q = o*(1 - v)*y + 26

    Y.append(y)
    A.append(a)
    K.append(k)
    I.append(i)
    N.append(n)
    E.append(e)
    Eland.append(eland)
    Mat.append(mat)
    Muo.append(muo)
    Mlo.append(mlo)
    Tat.append(tat)
    Tlo.append(tlo)
    O.append(o)
    V.append(v)
    Q.append(q)
    F.append(f)
    S.append(s)
    US.append(us)
    LA.append(la)
    NR.append(nr)

#output data
data = pd.DataFrame({'Emissions':E, 'Output':Y, 'Capital':K, "Temperature":Tat})
data.to_csv('DICE_S_R.csv')

#plotting output
plt.figure(1)
plt.plot(Y[0:500])
plt.title("Output over time")
plt.xlabel("Year")
plt.ylabel("Trillions of Dollars USD")
plt.savefig("DICE S_R_Y.png")

#plotting emissions
plt.figure(2)
plt.plot(E[0:500])
plt.title("Emissions over time")
plt.xlabel("Year")
plt.ylabel("Emissions CO2 in ppm")
plt.savefig("DICE S_R_E.png")

#plotting resourcess
plt.figure(3)
plt.plot(NR[0:500])
plt.title("Natural Resources over time")
plt.xlabel("Year")
plt.ylabel("Exajoules")
plt.savefig("DICE S_R_NR")

plt.figure(6)
plt.plot(S[0:500])
plt.title("Renewable Energy Reserves")
plt.ylabel("Exajoules")
plt.xlabel("Year")
plt.savefig("DICE_S_R_S")

plt.figure(4)
plt.plot(US[0:500])
plt.title("Non-Renewable Energy Reserves")
plt.ylabel("Exajoules")
plt.xlabel("Year")
plt.savefig("DICE_S_R_US")

plt.figure(5)
plt.plot(LA[0:500])
plt.title("Land Development over Time")
plt.ylabel("Hectares of Undeveloped Land")
plt.xlabel("Year")
plt.savefig("DICE_S_R_LA")

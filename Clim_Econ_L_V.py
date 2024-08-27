#Imports
import pandas as pd
import numpy as np

#Inital values
T=500#num years
GWP_0=67.79 #inital gwp
alpha=.005 #damage to gwp
beta=.03 #increase to gwp
eps=.1 #emissions/gwp
B_0=1527 #biosphere initial co2
Mlo_0=10010 #inital lower ocean c02
GHG_0=280 #historical value of ghg
GHG_1=830 #inital value of ghg
cs=3.2 #climate sensitivity
T_0=np.log(2)*np.log(GHG_1/GHG_0)*cs #initial avg temp

#Equations
#T(t)=ln(2)*ln(GHG(t)/GHG(0))*cs
#GHG(t+1)=.988*GHG(t)+.0047*(B(t))+eps*GWP(t)
#B(t+1)=.9948*(B(t))+.012*GHG(t)+.0001*Mlo(t)
#Mlo(t+1)=.9999*Mlo(t)+.0005*B(t)
#GWP(t+1)=(1+beta -alpha*T(GHG(t)))*GWP(t)

#Arrays to hold model values
GHG_vals=[]
B_vals=[]
Mlo_vals=[]
GWP_vals=[]
T_vals=[]

#calc model
for t in range(0,T):
    if(t==0):
        #use inital values to calc year 1
        ghg = 0.988*GHG_1 + 0.0047*B_0 + eps*GWP_0
        gwp = GWP_0 + GWP_0*beta - alpha*T_0*GWP_0
        temp = np.log(2)*np.log(GHG_1/GHG_0)*cs
        b = 0.9948*B_0 + 0.012*GHG_1 + 0.0001*Mlo_0
        mlo = 0.9999*Mlo_0 + 0.0005*B_0
    else:
        #use previous values to calculate current year
        ghg = 0.988*GHG_vals[t-1] + 0.0047*B_vals[t-1] + eps*GWP_vals[t-1]
        gwp = GWP_vals[t-1] + GWP_vals[t-1]*beta - alpha*T_vals[t-1]*GWP_vals[t-1]
        temp = np.log(2)*np.log(GHG_vals[t-1]/GHG_0)*cs
        b = 0.9948*(B_vals[t-1]) + 0.012*GHG_vals[t-1] + 0.0001*Mlo_vals[t-1]
        mlo = 0.9999*Mlo_vals[t-1] + 0.0005*B_vals[t-1]
  #Append values to respective arrays
    GHG_vals.append(ghg)
    B_vals.append(b)
    Mlo_vals.append(mlo)
    GWP_vals.append(gwp)
    T_vals.append(temp)

data = pd.DataFrame({'GHG':GHG_vals, 'B':B_vals, 'Mlo':Mlo_vals, 'GWP':GWP_vals,'delta-Temp':T_vals})

data.to_csv('Climate_Economy_L_V.csv')

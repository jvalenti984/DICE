#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 20:51:32 2023

@author: jarrettvalenti
"""
import numpy as np
import matplotlib.pyplot as plt

T = 200
E_d = []
E_c = []
E_0=830
M= 588+1350
Y= 67.78

for i in range(T):
    if i==0:
        e_d = E_0
        e_c = E_0*np.e**(-0.012*i) + 13095.4*np.e**(0.012*i)
    else:
        e_d = .988*E_d[i-1]+0.0047*M*i + 0.1*Y*i
        e_c = E_0*np.e**(-0.012*i) + 13095.4*np.e**(0.012*i)a
    E_d.append(e_d)
    E_c.append(e_c)

for e in range(len(E_c)):
    E_c[e]-=13095.4

plt.figure(1)
plt.plot(E_d)
plt.plot(E_c)
plt.title("Continuous vs Discrete Emissions Functions")
plt.xlabel("Year")
plt.ylabel("ppm CO2")
plt.legend(["Discrete Emissions", "Continuous Emissions"])
plt.savefig("cont_disc_emissions_200.png")

pct=[]
for i in range(T):
    pct.append(np.abs(E_d[i]-E_c[i])/(E_c[i]) * 100)

plt.figure(2)
plt.plot(pct)
plt.title("Percent Error between Emissions Functions")
plt.xlabel("Year")
plt.ylabel("Percent Error")
plt.savefig("pcte_disc_cont_emissions_200.png")

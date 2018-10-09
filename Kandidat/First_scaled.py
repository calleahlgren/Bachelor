#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 15:02:58 2018

@author: calle
"""

from math import *
import matplotlib.pyplot as plt


M = 0.01*3.003*10**(-6) #Earth masses * Scaling to Solar masses
r = 25 #AU 
r_vals = [r]
M_vals = [M]
#M_vals.append(M)
#r_vals.append(r)

beta = 1
zeta = 3/7
k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)
M_g= 10**-8 #Solar masses/yr
G = 4*pi**2 #AU**3/Solar Masses*yr**2 
M_sol = 1 #Solar masses
c_s0 = 2.04984*10**10/149597870700 #AU/yr
v_r = -31536000/149597870700  #AU/yr
u_r = v_r
St = 0.01
xi = 0.01
M_p = M_g*xi
delta_t = 100 # yr
t=0

Sigma_g = -M_g/(2*pi*r*u_r)
print(Sigma_g)

#Loop here

while r > 0.1:    

    Sigma_g = -M_g/(2*pi*r*u_r)

    Sigma_p = -M_p/(2*pi*r*v_r)

    Omega = sqrt(G*M_sol/r**3)

    R_H = ((M/(3*M_sol))**(1/3))*r

    v_k = sqrt(G*M_sol/r)

    c_s = c_s0*r**(-zeta/2)

    H = c_s/Omega

    M_delta = 2*((St/0.1)**(2/3))*Omega*R_H**2*Sigma_p

    r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(H/r)**(-2)*v_k

    M += M_delta*delta_t 

    r += r_delta*delta_t
    
    t += delta_t
    
    M_vals.append(M)
    r_vals.append(r)
    
for i in range(len(M_vals)):
    M_vals[i] *= 332946

plt.ylim(0.01, 1000)
plt.xlim(0.1, 100)
plt.loglog(r_vals,M_vals)
print(t)



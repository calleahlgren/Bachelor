#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 17:25:31 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

M_sol = 1.9885*10**30
G = 6.67408*10**(-11)
AU = 1.496*10**11
M_E = 5.972*10**24
r = 25*AU
M = 0.01*M_E
year = 3600*24*365

St = 0.01
alpha = 0.01
alpha_v = alpha*0.1
beta = 15/14
zeta = 3/7
gamma = (3/2) - zeta
chi = beta + (zeta/2) + (3/2)
xi = 0.02

k_mig = 2*(1.36 + (0.62*beta) + (0.43*zeta))

t = year*0.9*10**6
dt = 100*year

timeStamp = year*1*10**6

M_g_Ini = (M_sol*10**(-7))/year
M_g_Fin = (M_sol*10**(-8))/year

t_s = (year*3*10**6)/((M_g_Fin/M_g_Ini)**(-(2-gamma)/(5/2 - gamma)) - 1)

r_vals = [r]
M_vals = [M]
iso_vals = []

tR_vals = []
tM_vals = []

while t < year*3*10**6:
    
    omega = (G*M_sol/(r**3))**(1/2)
    v_k = ((G*M_sol)/r)**(1/2)
    R_H = r*((M/(3*M_sol))**(1/3))
    cS = 650*(r/AU)**(-zeta/2)
    H = (cS/omega)
    
    u_r = -(3/2)*alpha*cS*(H/r)
    v_delta = (1/2)*(H/r)*chi*cS
    v_r = (-2*St*v_delta) + u_r
    
    M_g = M_g_Ini*((t/t_s) + 1)**(-((5/2)-gamma)/(2-gamma))
    M_p = xi*M_g
    
    Sigma_g = -M_g/(2*pi*r*u_r)
    Sigma_p = -M_p/(2*pi*r*v_r)
    
    M_iso = 25*M_E*(1-((-chi+2.5)/6))*((0.34*((log(0.001, 10)/log(alpha_v, 10))**4)) + 0.66)*(((H/r)/0.05)**3)
    
#    if M < M_iso:
    M_delta = 2*((St/0.1)**(2/3))*omega*Sigma_p*(R_H**2)
#    a = 2*((St/0.1)**(2/3))
#    b = G*M_sol*(3*M_sol)**(-2/3)
#    c = 2*pi*((chi*St) + ((3/2)*alpha))*650**2*AU**(zeta)
#    d = xi*M_g/c
#    e = (M**(2/3))*(r**(zeta-1))
#    M_delta = a*b*d*e
        
    r_delta = -k_mig*(M/M_sol)*((Sigma_g*(r**2))/M_sol)*((H/r)**(-2))*v_k
    
    r_delta = r_delta/(1+(M/(2.3*M_iso))**2)
    
    M = M + (M_delta*dt)
    r = r + (r_delta*dt)
    
    r_vals.append(r)
    M_vals.append(M)
    iso_vals.append(M_iso/M_E)
    
    t = t+dt
    
    if t == timeStamp:
        tM_vals.append(M/M_E)
        tR_vals.append(r/AU)
        timeStamp += year*10**6
    

for i in range(len(M_vals)):
    M_vals[i] = M_vals[i]/M_E
    r_vals[i] = r_vals[i]/AU
    
    
plt.figure(figsize=(10,4))
plt.ylim(0.01, 3000)
plt.xlim(0.1, 100)
plt.loglog(r_vals,M_vals)
plt.loglog(tR_vals,tM_vals, 'ro')
plt.loglog(r_vals[1:],iso_vals)
plt.savefig('NewRun.png', dpi=500)

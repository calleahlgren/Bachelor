#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 02:33:27 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
from Bachelor_Iso_Mass import isoMass

r = np.linspace(25, 0.1, 1000)
r1 = np.linspace(5, 0.1, 1000)
r2 = np.linspace(100, 0.1, 1000)
M_e = 5.972*10**24

def analyticalGrowth(r,alpha,St,xi,M_0,r_0):
    
    M_0 = M_0*3.003*10**(-6) #Earth masses * Scaling to Solar masses
#    M_0 = 0.01*3.003*10**(-6) #Earth masses * Scaling to Solar masses
#    r_0 = 25 #AU
    
#    alpha = 0.01
    beta = 1
    zeta = 3/7
    chi = beta + (zeta/2) + (3/2) 
    k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)
    M_g= -10**-8 #Solar masses/yr
    G = 4*pi**2 #AU**3/Solar Masses*yr**2 
    M_sol = 1 #Solar masses
    c_s0 = 2.04984*10**10/149597870700 #AU/yr
    v_r = 31536000/149597870700  #AU/yr
    u_r = v_r
#    St = 0.01
#    xi = 0.01
    M_p = M_g*xi
    
    a = (((4/3)*xi)/((2/3)*(St/alpha)*chi + 1))
    b = ((2*(St/0.1)**(2/3)*M_sol*(M_sol)**(-2/3))/(k_mig*G*c_s0**(-2)))
    c = (r**(1-zeta) - r_0**(1-zeta))/(1-zeta)
    
    M = M_0**(4/3) - (a * b * c)
    M = M**(3/4)
    
    return M*332946

plt.ylim(0.01, 1000)
plt.xlim(0.1, 100)
plt.loglog(r,M(r, 0.01, 0.01, 0.01, 0.01, 25))
plt.loglog(r,M(r, 0.01, 0.01, 0.02, 0.01, 25))
plt.loglog(r,M(r, 0.001, 0.001, 0.01, 0.01, 25))
plt.loglog(r,M(r, 0.001, 0.001, 0.02, 0.01, 25))
#plt.loglog(r1,M(r1, 0.01, 0.01, 0.01, 0.01, 5))
#plt.loglog(r1,M(r1, 0.01, 0.01, 0.02, 0.01, 5))
#plt.loglog(r1,M(r1, 0.001, 0.001, 0.01, 0.01, 5))
#plt.loglog(r1,M(r1, 0.001, 0.001, 0.02, 0.01, 5))
plt.loglog(r2,isoMass(r2, 0.01)/M_e)
plt.loglog(r2,isoMass(r2, 0.001)/M_e)
print(M(r[-1], 0.01, 0.01, 0.01, 0.01, 25))

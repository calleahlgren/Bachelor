#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 03:29:05 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

def isoMass(r, alpha):
    
    M_e = 5.972*10**24
    AU = 149597870700
    alpha_v = 0.1*alpha    
    beta = 1
    zeta = 3/7
    chi = beta + (zeta/2) + (3/2)
    c_s0 = 650
    G = 6.67408*10**-11  
    M_sol = 1.989*10**30
    r = r*AU
    
    c_s = c_s0*(r/AU)**(-zeta/2)
    
#    Omega = sqrt(G*M_sol/r**3)
#    Omega = G*M_sol/(r**3)
    Omega = (G*M_sol/r**3)**(1/2)
    
    H = c_s/Omega    
    
    a = 25*M_e*((H/r)/0.05)**3
#    a = 25*((H/r)/0.05)**3
    b = (log(alpha, 10)/log(alpha_v, 10))**4
#    b = (log(alpha)/log(alpha_v))**4
#    b = (log(0.001, 10)/log(alpha_v, 10))**4
    c = (0.34*b + 0.66)
    d = (1+((chi + 2.5)/6))
    
    M_iso = a*c*d
#    M_iso /= M_e
    return M_iso

#M_iso_0 = isoMass(25, 0.01)
#print(M_iso_0)
#AU = 149597870700
#r = np.linspace(25*AU, 0.1*AU, 25000)
#r_scaled = [i/AU for i in r]
#plt.loglog(r_scaled,isoMass(r, 0.01))
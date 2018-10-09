#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 17:58:20 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

def growthTrack(St, alpha, xi, M_0, r_0):

    AU = 149597870700
    year = 24*3600*365
    M_e = 5.972*10**24
    M = M_0*M_e
    r = r_0*AU # 
    r_vals = [r]
    M_vals = [M]

    beta = 1
    zeta = 3/7
    k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)
    M_g= (10**-8)*(1.989*10**30)/year
    G = 6.67408*10**-11  
    M_sol = 1*1.989*10**30
    c_s0 = 650
    M_p = M_g*xi
    delta_t = year*10
    t=0
    
    while r > 10000:    
    
        Omega = sqrt(G*M_sol/r**3)
    
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = sqrt(G*M_sol/r)
    
        c_s = c_s0*(r/AU)**(-zeta/2)
    
        H = c_s/Omega
        
        u_r = -(3/2)*alpha*c_s*(H/r)
        v_r = u_r
        
        Sigma_g = -M_g/(2*pi*r*u_r)
    
        Sigma_p = -M_p/(2*pi*r*v_r)
    
        M_delta = 2*((St/0.1)**(2/3))*Omega*R_H**2*Sigma_p
    
        r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(H/r)**(-2)*v_k
    
        M += M_delta*delta_t 
    
        r += r_delta*delta_t
        
        t += delta_t
        
        M_vals.append(M)
        r_vals.append(r)
        
    for i in range(len(M_vals)):
        M_vals[i] /= M_e
        
    for i in range(len(r_vals)):
        r_vals[i] /= AU
        
    return r_vals, M_vals 

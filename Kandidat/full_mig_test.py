#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 00:14:06 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

M_e = 5.972*10**24
AU = 149597870700
beta = 1
zeta = 3/7
chi = beta + (zeta/2) + (3/2)
c_s0 = 650
G = 6.67408*10**-11  
M_sol = 1.989*10**30
k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)
year = 24*3600*36
G = 6.67408*10**-11

def omega(r):
    
    return (G*M_sol/r**3)**(1/2)

def cS(r):
    
    return c_s0*(r/AU)**(-zeta/2)

def scaleHeight(r):
    
    return cS(r)/omega(r)

def isoMass(r, alpha): 
    
    alpha_v = 0.1*alpha

#    r *= AU
#    r_1 = r*AU
    r = r*AU
    
    a = 25*M_e*((scaleHeight(r)/r)/0.05)**3
    b = (log(alpha, 10)/log(alpha_v, 10))**4
    c = (0.34*b + 0.66)
    d = (1+((chi + 2.5)/6))
    
    M_iso = a*c*d
#    M_iso /= M_e
    return M_iso

#r = np.linspace(100, 0.1, 25000)
#print(r[0])
#
#plt.ylim(0.01, 1000)
#plt.xlim(0.1, 100)
#print(r[0])
#
#plt.loglog(r,(isoMass(r, 0.01)/M_e))
#print(r[0])

def maxMass(St, xi, alpha, r_0):
    
    r_0 = r_0*AU
    
    a = (2.9*(St/0.01)**(2/3)/((2/3)*(St/alpha)*chi+1))**(3/4)
    b = (xi/0.01)**(3/4)
    c = (k_mig/4.42)**(-3/4)
    d = ((1-zeta)/(4/7))**-1
    e = (r_0/(25*AU))**((3/4)*(1-zeta))
    
    M_max = 11.7*M_e*a*b*c*d*e
    
    return M_max

def growthTrack(alpha, St, xi, M_0, M_g, r_0):    
    
    r = r_0*AU
    M = M_0*M_e    
    
    v_r = -1
    u_r = v_r
    M_g = M_g*M_sol/year
    M_p = M_g*xi
    
    r_vals = [r]
    M_vals = [M]
    
    delta_t = year 
    t=0
    
    M_iso_0 = isoMass(r_0, alpha)
    print(M_iso_0)
    
    def mKH(M):
    
     kappa = 0.005
     
     M_kh = ((M_e*10**-5)/year)*(M/(10*M_e))**4*(kappa/0.1)**-1
     
     return M_kh 
    
    def mDisc(r, alpha, M_g, M):
        
        a = (1.5*10**-3*M_e/year)
        b = ((scaleHeight(r)/r)/0.05)**-4
        c = (M/(10*M_e))**(4/3)
        d = (alpha/0.01)**-1
        e = (M_g/((10**-8)*M_sol/year))
        f = (1/(1+(M/(M_iso_0/1.5))**2))
        
        M_disc = a*b*c*d*e*f
        
        return M_disc
    
    def mG(r, alpha, M_g, M):
        
        return min(mKH(M), mDisc(r, alpha, M_g, M), M_g)

    while t < year*5.5*10**6:
        
        Sigma_g = -M_g/(2*pi*r*u_r)
        
        Sigma_p = -M_p/(2*pi*r*v_r)
    
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)
        
#        M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
        
        if M < M_iso_0:   
            
            M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
        
        else:
            
            M_delta = mG(r, alpha, M_g, M)
#            
#        print(M_delta)
        
        r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(scaleHeight(r)/r)**(-2)*v_k
            
        r_delta /= (1 + (M/(1.5*isoMass((r/AU), alpha)))**2)
            
        M += M_delta*delta_t
        
#        print(M)
        
        r += r_delta*delta_t    
        
        t += delta_t
            
        M_vals.append(M)
        r_vals.append(r)
            
    for i in range(len(M_vals)):
        M_vals[i] /= 5.972*10**24
    
    for i in range(len(r_vals)):
        r_vals[i] /= AU
    
#    plt.ylim(0.01, 1000)
#    plt.xlim(0.1, 100)
    plt.loglog(r_vals, M_vals)
#    print(t/year)
    
#    return plt.loglog(r_vals,M_vals, label='(alpha, St, xi)=(%.3f, %.3f, %.3f)' %(alpha, St, xi))
    return
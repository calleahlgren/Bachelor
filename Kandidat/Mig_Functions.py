#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 20:20:30 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

from Helper_Funcs import *

def maxMass(St, xi, alpha, r_0):
    
    r_0 = r_0*AU
    
    a = (2.9*(St/0.01)**(2/3)/((2/3)*(St/alpha)*chi+1))**(3/4)
    b = (xi/0.01)**(3/4)
    c = (k_mig/4.42)**(-3/4)
    d = ((1-zeta)/(4/7))**-1
    e = (r_0/(25*AU))**((3/4)*(1-zeta))
    
    M_max = 11.7*M_e*a*b*c*d*e
    
    return M_max


def analyticalGrowth(r,alpha,St,xi,M_0,r_0):
    
    AU = 149597870700
    
    M_0 = M_0*M_e 
    
    r = r*AU
    r_0 = r_0*AU
    
    a = (((4/3)*xi)/((2/3)*(St/alpha)*chi + 1))
    b = ((2*(St/0.1)**(2/3)*M_sol*(M_sol)**(-2/3))/(k_mig*G*c_s0**(-2)*AU**(-zeta)))
    c = (r**(1-zeta) - r_0**(1-zeta))/(1-zeta)
    
    M = M_0**(4/3) - (a * b * c)
    M = M**(3/4)
    
    return M

def growthTrack(St, alpha, xi, M_0, r_0):

    M = M_0*M_e
    r = r_0*AU # 
    r_vals = [r]
    M_vals = [M]

    M_g= (10**-8)*(1.989*10**30)/year
    M_p = M_g*xi
    delta_t = year*10 
    t=0
    
    while r > 10000:    
    
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)

        u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
        v_r = u_r
        
        Sigma_g = -M_g/(2*pi*r*u_r)
    
        Sigma_p = -M_p/(2*pi*r*v_r)
    
        M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
    
        r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(scaleHeight(r)/r)**(-2)*v_k
    
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

def growthTrack2(St, alpha, xi, M_0, r_0):

    M = M_0*M_e
    r = r_0*AU # 
    r_vals = [r]
    M_vals = [M]

    M_g= (10**-8)*(1.989*10**30)/year
    M_p = M_g*xi
    delta_t = year
    t=0
    
    while t < year*2.1*10**6:    
        
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)

        u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
        v_r = u_r
        
        Sigma_g = -M_g/(2*pi*r*u_r)
    
        Sigma_p = -M_p/(2*pi*r*v_r)
    
        M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
    
        r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(scaleHeight(r)/r)**(-2)*v_k
        
        r_delta /= (1+(M/(2.3*isoMass(r/AU, alpha)))**2)
    
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

def growthTrack3(St, alpha, xi, M_0, r_0):

    M = M_0*M_e
    r = r_0*AU # 
    r_vals = [r]
    M_vals = [M]

#    M_g= (10**-8)*(1.989*10**30)/year
#    M_p = M_g*xi
    delta_t = 5000*year
    t=year*0.9*10**6
    
    
#    M_iso_0 = isoMass(r_0, alpha)
#    M_max = maxMass(St, xi, alpha, r_0)
#    A = (M_iso_0/M_max)**(4/3)
#    X = ((1+4*A)**(1/2)-1)/(2*A)
#    r_iso = X**(1/(1-zeta))*r_0
#    
#    M_iso = M_iso_0*(r_iso/r_0)**(2*(1-zeta)*3/4)
#    M_iso = isoMass(r_iso, alpha)
#    
    
    """I used this to test out the timestep you mentioned, the next comment 
    will explain why"""
#    while t < year*3*10**6 and r>1*AU:    
    while t < year*3*10**6:
        
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)

        u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
#        v_r = u_r
        delta_v_r = 1/2*(scaleHeight(r)/r)*chi*cS(r)
        v_r = -2*St*delta_v_r + u_r
#        v_r = -(2*delta_v_r)/(St + St**-1) + u_r/(1 + St**2)
        
        M_g = M_g_Ini*(t/t_s + 1)**-((5/2-gamma)/(2-gamma))
        
        M_p = M_g*xi
        
        Sigma_g = -M_g/(2*pi*r*u_r)
    
        Sigma_p = -M_p/(2*pi*r*v_r)
    

        if M < isoMass(r/AU, alpha):
            
            M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
        else:
            
            M_delta = mG(r, alpha, M_g, M)
            
        
        r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(scaleHeight(r)/r)**(-2)*v_k
        
        r_delta /= (1+(M/(1.5*isoMass(r/AU, alpha)))**2)
#        r_delta /= (1+(M/(2.3*isoMass(r/AU, alpha)))**2)
        
        """I tried implementing this as the timestep but that made the steps
        far to big and all planets regardless of St, alpha and xi fell in to 
        the sun"""
#        delta_t = min(M/M_delta, abs(r/r_delta))
   
        M += M_delta*delta_t 
    
        r += r_delta*delta_t
        
        t += delta_t
#        
        if delta_t > 500*year:
        
            delta_t -= 100*year
#        
        M_vals.append(M)
        r_vals.append(r)
        
        if r<0.0001*AU:
            break
        
    for i in range(len(M_vals)):
        M_vals[i] /= M_e
        
    for i in range(len(r_vals)):
        r_vals[i] /= AU
        
    return r_vals, M_vals

    
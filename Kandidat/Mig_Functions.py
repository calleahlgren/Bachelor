#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 20:20:30 2018

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

def sigmaG(M_g, alpha, r):
    
    r = r*AU
    M_g= M_g*M_sol/year
    
    u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
    
    return -M_g/(2*pi*r*u_r)

def isoMass(r, alpha): 
    
    alpha_v = 0.1*alpha

    r = r*AU
    
    a = 25*M_e*((scaleHeight(r)/r)/0.05)**3
    b = (log(alpha, 10)/log(alpha_v, 10))**4
    c = (0.34*b + 0.66)
    d = (1+((chi + 2.5)/6))
    
    M_iso = a*c*d
    
    return M_iso

def maxMass(St, xi, alpha, r_0):
    
    r_0 = r_0*AU
    
    a = (2.9*(St/0.01)**(2/3)/((2/3)*(St/alpha)*chi+1))**(3/4)
    b = (xi/0.01)**(3/4)
    c = (k_mig/4.42)**(-3/4)
    d = ((1-zeta)/(4/7))**-1
    e = (r_0/(25*AU))**((3/4)*(1-zeta))
    
    M_max = 11.7*M_e*a*b*c*d*e
    
    return M_max

def mKH(M):
    
     kappa = 0.005
     
     M_kh = ((M_e*10**-5)/year)*(M/(10*M_e))**4*(kappa/0.1)**-1
     
     return M_kh 
    
def mDisc(r, alpha, M_g, M):
    
#    a = (1.5*10**-3*M_e/year)
#    b = ((scaleHeight(r)/r)/0.05)**-4
#    c = (M/(10*M_e))**(4/3)
#    d = (alpha/0.01)**-1
#    e = (M_g/((10**-8)*M_sol/year))
#    f = (1/(1+(M/(M_iso_0/1.5))**2))
#    
#    M_disc = a*b*c*d*e*f
    
    alpha_v = alpha*0.1
#    
    a = 0.29/(3*pi)
    b = (scaleHeight(r)/r)**(-4)
    c = (M/M_sol)**(4/3)
    d = M_g/alpha
    e = (M/M_sol)**2*(scaleHeight(r)/r)**(-5)*alpha_v**(-1)
    f = 1/(1+0.04*e)
    
    M_disc = a*b*c*d*f
    
    return M_disc
    
def mG(r, alpha, M_g, M):
    
    return min(mKH(M), mDisc(r, alpha, M_g, M), M_g)

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
        
        r_delta /= (1+(M/(1.5*isoMass(r/AU, alpha)))**2)
    
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

    M_g= (10**-8)*(1.989*10**30)/year
    M_p = M_g*xi
    delta_t = year
    t=0
    
#    M_iso_0 = isoMass(r_0, alpha)
#    M_max = maxMass(St, xi, alpha, r_0)
#    A = (M_iso_0/M_max)**(4/3)
#    X = ((1+4*A)**(1/2)-1)/(2*A)
#    r_iso = X**(1/(1-zeta))*r_0
#    
#    M_iso = M_iso_0*(r_iso/r_0)**(2*(1-zeta)*3/4)
#    M_iso = isoMass(r_iso, alpha)
#    
    
    
    while t < year*2.1*10**6:    
        
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)

        u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
        v_r = u_r
        
        Sigma_g = -M_g/(2*pi*r*u_r)
    
        Sigma_p = -M_p/(2*pi*r*v_r)
    
#        if M < M_iso:
        if M < isoMass(r/AU, alpha):
            
            M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
        
        else:
            
            M_delta = mG(r, alpha, M_g, M)
    
        r_delta = -k_mig*(M/M_sol**2)*Sigma_g*r**2*(scaleHeight(r)/r)**(-2)*v_k
        
        r_delta /= (1+(M/(1.5*isoMass(r/AU, alpha)))**2)
    
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


    
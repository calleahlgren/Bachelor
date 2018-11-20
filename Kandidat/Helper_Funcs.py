#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 18:47:35 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt


"""Constants and helper functions I use for both the growth track generating
function as well as the growth map generating one. All constants and variables
are named like they are named in your paper so everything should be reasonable
easy to understand. You'll find that in some of the functions I have broken
up the equations in to smaller parts, this is just to make it more readable and
easy to keep track of for myself"""

M_e = 5.972*10**24 # kg
AU = 149597870700 # m
beta = 15/14
zeta = 3/7
chi = beta + (zeta/2) + (3/2)
gamma = 15/14
c_s0 = 650 # m/s
G = 6.67408*10**-11 #  
M_sol = 1.989*10**30 # kg
k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)
year = 24*3600*365 # s

M_g_Ini = 10**-7*M_sol/year # kg/s
M_g_Fin = 10**-8*M_sol/year # kg/s

t_1 = 3*10**6*year # s
t_s = t_1/((M_g_Fin/M_g_Ini)**-((2-gamma)/(5/2-gamma)) - 1)

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

'Equation 24'
def isoMass(r, alpha): 
    
    alpha_v = 0.1*alpha

    r = r*AU
    
    a = 25*M_e*((scaleHeight(r)/r)/0.05)**3
    b = (log(0.001, 10)/log(alpha_v, 10))**4
    c = (0.34*b + 0.66)
    d = (1+((chi - 2.5)/6))
    
    M_iso = a*c*d
    
    return M_iso

'Equation 38'
def mKH(M):
    
     kappa = 0.005
     
     M_kh = ((M_e*10**-5)/year)*(M/(10*M_e))**4*(kappa/0.1)**-1
     
     return M_kh 
    
'Equation 39'    
def mDisc(r, alpha, M_g, M):
#    
    a = (1.5*10**-3*M_e/year)
    b = ((scaleHeight(r)/r)/0.05)**-4
    c = (M/(10*M_e))**(4/3)
    d = (alpha/0.01)**-1
    e = (M_g/((10**-8)*M_sol/year))
#    f = (1/(1+(M/(1.5*isoMass(r/AU, alpha)))**2))
    f = (1/(1+(M/(2.3*isoMass(r/AU, alpha)))**2))
    
    M_disc = a*b*c*d*e*f
    """After I fixed the problem with the isolation mass function the difference
    between these two version of Eq 39 became a lot less apparent, but there
    is still a slight difference."""
#    alpha_v = alpha*0.1
##    
#    a = 0.29/(3*pi)
#    b = (scaleHeight(r)/r)**(-4)
#    c = (M/M_sol)**(4/3)
#    d = M_g/alpha
#    e = (M/M_sol)**2*(scaleHeight(r)/r)**(-5)*alpha_v**(-1)
#    f = 1/(1+0.04*e)
#    
#    M_disc = a*b*c*d*f
    
    return M_disc
    
'Equation 40'
def mG(r, alpha, M_g, M):
    
    return min(mKH(M), mDisc(r, alpha, M_g, M), M_g)
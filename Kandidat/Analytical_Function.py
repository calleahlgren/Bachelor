#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 02:33:27 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

def analyticalGrowth(r,alpha,St,xi,M_0,r_0):
    
    AU = 149597870700
    M_e = 5.972*10**24
    
    M_0 = M_0*M_e 
    
    r = r*AU
    r_0 = r_0*AU
    
    beta = 1
    zeta = 3/7
    chi = beta + (zeta/2) + (3/2) 
    k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)    
    G = 6.67408*10**-11  
    M_sol = 1*1.989*10**30
    c_s0 = 650
    
    a = (((4/3)*xi)/((2/3)*(St/alpha)*chi + 1))
    b = ((2*(St/0.1)**(2/3)*M_sol*(M_sol)**(-2/3))/(k_mig*G*c_s0**(-2)*AU**(-zeta)))
    c = (r**(1-zeta) - r_0**(1-zeta))/(1-zeta)
    
    M = M_0**(4/3) - (a * b * c)
    M = M**(3/4)
    
    return M


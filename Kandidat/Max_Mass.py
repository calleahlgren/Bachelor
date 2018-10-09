#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 02:45:34 2018

@author: calle
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt


def maxMass(St, xi, alpha, r_0):
    
    AU = 149597870700
    M_e = 5.972*10**24
    M_sol = 1.989*10**30
    beta = 1
    zeta = 3/7
    chi = beta + (zeta/2) + (3/2)
    c_s0 = 650
    k_mig = 2*(1.36 + 0.62*beta + 0.43*zeta)
    r_0 *= AU
    
    a = (2.9*(St/0.01)**(2/3)/((2/3)*(St/alpha)*chi+1))**(3/4)
    b = (xi/0.01)**(3/4)
    c = (k_mig/4.42)**(-3/4)
    d = ((1-zeta)/(4/7))**-1
    e = (r_0/(25*AU))**((3/4)*(1-zeta))
    
    M_max = 11.7*M_e*a*b*c*d*e
    
    return M_max
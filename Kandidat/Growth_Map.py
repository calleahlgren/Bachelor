#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 18:42:59 2018

@author: calle
"""

from Helper_Funcs import*

"""This function is identical to the growth track generating function aside
from the output which is in format that is better for generating growth maps."""

def growthMap(t_0, St, alpha, xi, M_0, r_0):

    M = M_0*M_e
    r = r_0*AU # 
    
    'Time step'
    delta_t = 100*year
    
    'Starting time'
    t = t_0*year*10**6
    
    
    while t < year*3*10**6:    
        
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)

        u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
        delta_v_r = 1/2*(scaleHeight(r)/r)*chi*cS(r)
        v_r = -2*St*delta_v_r + u_r

        
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
    
        M += M_delta*delta_t 
    
        r += r_delta*delta_t
        
        t += delta_t
        
#        if delta_t > 500*year:
#        
#            delta_t -= 100*year
        
        if r<0.001*AU:
            break
        
    return (log(M/M_e, 10), r/AU)
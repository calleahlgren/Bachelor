#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 21:15:55 2018

@author: calle
"""

from Helper_Funcs import *

"""This is the function I use to generate the growth tracks. I have tried to 
name all variables as they are named in the paper so I think the code should
be quite self explanatory"""

def growthTrack(St, alpha, xi, M_0, r_0):

    M = M_0*M_e
    r = r_0*AU # 
    r_vals = [r]
    M_vals = [M]
    
    timeStampsTime1 = 1.2*10**6*year
    timeStampsTime2 = year*10**6
    xTimeStamps1 = []
    yTimeStamps1 = []
    xTimeStamps2 = []
    yTimeStamps2 = []
    
    'Time step'
#    delta_t = 1000*year
    
    'Starting time'
    t=year*0.9*10**6
    
    
    while t < year*3*10**6:
        
        R_H = ((M/(3*M_sol))**(1/3))*r
    
        v_k = (G*M_sol/r)**(1/2)

        u_r = -(3/2)*alpha*cS(r)*(scaleHeight(r)/r)
        delta_v_r = 1/2*(scaleHeight(r)/r)*chi*cS(r)
        v_r = -2*St*delta_v_r + u_r
#        v_r = -(2*delta_v_r)/(St + St**-1) + u_r/(1 + St**2)
        
        M_g = M_g_Ini*(t/t_s + 1)**-((5/2-gamma)/(2-gamma))
#        M_g = 10**-8*M_sol/year
        M_p = M_g*xi
        
        Sigma_g = -M_g/(2*pi*r*u_r)   
        Sigma_p = -M_p/(2*pi*r*v_r)
    
        if M < isoMass(r/AU, alpha):
            
            M_delta = 2*((St/0.1)**(2/3))*omega(r)*R_H**2*Sigma_p
        else:
            
            M_delta = mG(r, alpha, M_g, M)
            print(t/year)
            break
                    
        r_delta = -k_mig*(M/(M_sol**2))*Sigma_g*r**2*(scaleHeight(r)/r)**(-2)*v_k        
#        r_delta /= (1+(M/(1.5*isoMass(r/AU, alpha)))**2)
        """Using M_T = 2.3*M_iso as below, as opposed to M_T = 1.5*M_iso as 
        above, made the growth tracks further diverge from yours. Don't forget 
        to change this in the function called mDisc in Helper_Funcs as well if 
        you want to generate tracks with the 2.3 version"""
        r_delta /= (1+(M/(2.3*isoMass(r/AU, alpha)))**2)
        
        """I tried implementing this as the timestep but that made the steps
        far to big and all planets starting at 5 and 25 AU fell in to the Sun
        regardless of St, alpha and xi"""
        delta_t = 0.01*min(M/(M_delta), abs(r/(r_delta)))
   
        M += M_delta*delta_t 
    
        r += r_delta*delta_t
        
        t += delta_t
        
        if t > year*3.0*10**6:
            break
#        
#        if delta_t > 500*year:
#        
#            delta_t -= 100*year
#        
        M_vals.append(M)
        r_vals.append(r)
        
#        if t == 1.2*10**6*year:
#            print(r/AU)
#            print(t)
#        print(t)
#        if t == timeStampsTime1:
#            timeStampsTime1 += year*0.2*10**6
#            xTimeStamps1.append(r/AU)
#            yTimeStamps1.append(M/M_e)
#            print(t)
        if t == timeStampsTime2:
            timeStampsTime2 += year*10**6
            timeStampsTime1 += year*0.2*10**6
            xTimeStamps2.append(r/AU)
            yTimeStamps2.append(M/M_e)
#        
        if r<0.0001*AU:
            break
        
    for i in range(len(M_vals)):
        M_vals[i] /= M_e
        
    for i in range(len(r_vals)):
        r_vals[i] /= AU
        
#    print(M_vals[-1])
#    print(t)
        
    return r_vals, M_vals, xTimeStamps1, yTimeStamps1, xTimeStamps2, yTimeStamps2
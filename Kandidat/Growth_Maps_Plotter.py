#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:40:25 2018

@author: calle
"""

from Growth_Map import *
from matplotlib import ticker

"""Similar to Growth_Track_Plotter this is simply the code I put together to 
plot growth maps."""

r = np.linspace(0.1, 50, 100)

t = np.linspace(0, 3, 100)

X, Y = np.meshgrid(r, t)

xi_vals = [0.01, 0.02, 0.05, 0.1]

St = 1


for k in range(2):
    
    St /= 10
    alpha = St

    for l in range(4):
        
        Z1 = np.zeros((100,100))
        Z2 = np.zeros((100,100))
        
        
        for i in range(100):
            for j in range(100):
                cacheTuple = growthMap(Y[i,j],St, alpha, xi_vals[l], 0.01, X[i,j])
#                cacheTuple = growthMap(Y[i,j],0.01, 0.01, 0.01, 0.01, X[i,j])
                Z1[i,j] = cacheTuple[0]
                Z2[i,j] = cacheTuple[1]
        
        CS1 = plt.contourf(X, Y, Z1, 100, cmap=plt.get_cmap('gist_rainbow_r'))
#        plt.contour(X, Y, Z2, [1, 5, 10, 15, 20])
        #plt.annotate(r'$\xi = 0.01$', (0.78,0.95), xycoords = 'axes fraction')
        #plt.annotate(r'$\alpha$' f" = St = {St}0.0", (0.73,0.90), xycoords = 'axes fraction')
        plt.annotate(r'$\xi$' f" = {xi_vals[l]}" ,(0.78,0.95), xycoords = 'axes fraction')
        plt.annotate(r'$\alpha$' f" = St = {St}",(0.73,0.90), xycoords = 'axes fraction')
        
        csb = plt.colorbar(CS1, ticks=[-2, -1, 0, 1, 2, 3, 4])
        
        CS2 = plt.contour(X, Y, Z2, [1, 5.2, 10, 15, 20], colors = 'red')
#        CS3 = plt.contour(X, Y, Z1, [2.5], colors = 'black')

        plt.clabel(CS2, inline=1, fontsize=10)
#        plt.clabel(CS3, inline=1, fontsize=10, fmt = r'318 $M_E$')
        
        csb.set_label(r'log($M/M_E$)')
        plt.ylabel(r'$t_0$[Myr]')
        plt.xlabel(r'$r_0$[AU]')
#        plt.savefig('Growth_map1.png', dpi=500)
        
        if k == 0:
            plt.savefig('Growth_map%d.png' %(l+1), dpi=500)
            np.save(f'mass_data{l+1}', Z1)
            np.save(f'radial_data{l+1}', Z2)
        else:
            plt.savefig('Growth_map%d.png' %(l+5), dpi=500)
            np.save(f'mass_data{l+5}', Z1)
            np.save(f'radial_data{l+5}', Z2)
            
        plt.figure()

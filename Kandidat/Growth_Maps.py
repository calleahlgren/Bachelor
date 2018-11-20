#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:40:25 2018

@author: calle
"""
from Mig_Functions import *
from tempfile import TemporaryFile
from matplotlib import ticker

r = np.linspace(0.1, 50, 100)

t = np.linspace(0, 3, 100)

X, Y = np.meshgrid(r, t)

Z1 = np.zeros((100,100))
#Z2 = Z1[:]
Z2 = np.zeros((100,100))


for i in range(100):
    for j in range(100):
        cacheTuple = growthTrack4(Y[i,j],0.01, 0.01, 0.01, 0.01, X[i,j])
        Z1[i,j] = cacheTuple[0]
        Z2[i,j] = cacheTuple[1]

#for i in range(75):
#    for j in range(75):
#        Z2[i,j] = growthTrack4(Y[i,j],0.01, 0.01, 0.01, 0.01, X[i,j])[1]
        
cs = plt.contourf(X, Y, Z1)
plt.contour(X, Y, Z2, [1, 5, 10, 15, 20])
plt.annotate(r'$\Xi = \alpha$ = St = 0.01',(40,2.6))

csb = plt.colorbar(cs)
csb.set_label(r'log($M/M_E$)')
plt.ylabel(r'$t_0$[Myr]')
plt.xlabel(r'$r_0$[AU]')
plt.savefig('Growth_map2.png', dpi=500)


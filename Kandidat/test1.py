#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 16:20:38 2018

@author: calle
"""

import numpy as np
import matplotlib.pyplot as plt
from shapely import geometry


r = np.linspace(0.1, 50, 100)

t = np.linspace(0, 3, 100)

X, Y = np.meshgrid(r, t)

Z1 = np.load('mass_data5.npy')
Z2 = np.load('radial_data5.npy')

#massList = []
#iMList = []
#jMList = []
#
#radialList = []
#iRList = []
#jRList = []
#
#
#for i in range(100):
#    for j in range(100):
#        if x[i,j]>2.476 and x[i,j]<2.478:
#            iMList.append(i)
#            jMList.append(j)
#            massList.append(x[i,j])
#        if y[i,j]<5.22 and y[i,j]>5.18:
#            iRList.append(i)
#            jRList.append(j)
#            radialList.append(y[i,j])
#            
#   
#print(iMList)
#print(jMList)                 
#print(massList)
#print(iRList)
#print(jRList)                 
#print(radialList)


plt.contourf(X, Y, Z1, 100, cmap=plt.get_cmap('gist_rainbow_r'))
#CS2 = plt.contour(X, Y, Z2, [1, 5.2, 10, 15, 20], colors = 'red')
CS2 = plt.contour(X, Y, Z2, [5.2], colors = 'red')
CS3 = plt.contour(X, Y, Z1, [2.5], colors = 'black')
#fmt = {}
#strs = ['1 AU', '5.2 AU', '10 AU', '15 AU', '20 AU']
#for l, s in zip(CS2.levels, strs):
#    fmt[l] = s

#plt.clabel(CS2, inline=1, fontsize=10, fmt = fmt, manual = True)
#plt.clabel(CS2, inline=1, fontsize=10)
#plt.clabel(CS3, inline=1, fontsize=10, fmt = r'318 $M_E$')

def findIntersection(contour1,contour2):
  p1 = contour1.collections[0].get_paths()[0]
  v1 = p1.vertices

  p2 = contour2.collections[0].get_paths()[0]
  v2 = p2.vertices

  poly1 = geometry.LineString(v1)
  poly2 = geometry.LineString(v2)

  intersection = poly1.intersection(poly2)

  return intersection



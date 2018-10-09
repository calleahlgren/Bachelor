#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 21:18:56 2018

@author: calle
"""

from Mig_Functions import *

r = np.linspace(25, 0.1, 1000)
r1 = np.linspace(5, 0.1, 1000)
r2 = np.linspace(100, 0.1, 1000)

plt.figure(figsize=(10,4))
plt.ylim(0.01, 1000)
plt.xlim(0.1, 100)
gLine1, = plt.loglog(r,analyticalGrowth(r, 0.01, 0.01, 0.01, 0.01, 25)/M_e, 'b', label=r'($\rm \alpha$, St, $\rm \xi$)=0.01,0.01,0.01')
gLine2, = plt.loglog(r,analyticalGrowth(r, 0.01, 0.01, 0.02, 0.01, 25)/M_e, 'g', label=r'($\rm \alpha$, St, $\rm \xi$)=0.01,0.01,0.02')
gLine3, = plt.loglog(r,analyticalGrowth(r, 0.001, 0.001, 0.01, 0.01, 25)/M_e, 'r', label=r'($\rm \alpha$, St, $\rm \xi$)=0.001,0.001,0.01')
gLine4, = plt.loglog(r,analyticalGrowth(r, 0.001, 0.001, 0.02, 0.01, 25)/M_e, 'y', label=r'($\rm \alpha$, St, $\rm \xi$)=0.001,0.001,0.02')
plt.loglog(r1,analyticalGrowth(r1, 0.01, 0.01, 0.01, 0.01, 5)/M_e, 'b')
plt.loglog(r1,analyticalGrowth(r1, 0.01, 0.01, 0.02, 0.01, 5)/M_e, 'g')
plt.loglog(r1,analyticalGrowth(r1, 0.001, 0.001, 0.01, 0.01, 5)/M_e, 'r')
plt.loglog(r1,analyticalGrowth(r1, 0.001, 0.001, 0.02, 0.01, 5)/M_e, 'y')
iLine1, = plt.loglog(r2,isoMass(r2, 0.01)/M_e, '--', label=r'PIM $\rm \alpha_v$ = 0.001')
iLine2, = plt.loglog(r2,isoMass(r2, 0.001)/M_e, '--', label=r'PIM $\rm \alpha_v$ = 0.0001')

plt.ylabel(r'$M[M_{\rm E}]$')
plt.xlabel(r'$r$[AU]')
firstLegend = plt.legend(handles=[gLine1, gLine2, gLine3, gLine4], loc=3, fontsize='small', frameon=False, labelspacing=0.2)
ax = plt.gca().add_artist(firstLegend)
plt.legend(handles=[iLine1, iLine2], fontsize='small', loc=1, frameon=False, labelspacing=0.2)
plt.savefig('Analytical_growthTracks.png', dpi=500)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 01:35:18 2018

@author: calle
"""

from Growth_Track import *

"""This is simply the code I put together to plot growth tracks."""

r = np.linspace(100, 0.1, 25000)

plt.figure(figsize=(10,4))
plt.ylim(0.01, 3000)
plt.xlim(0.1, 100)
iLine1, = plt.loglog(r,(isoMass(r, 0.01)/M_e), '--', label=r'PIM $\rm \alpha_v$ = 0.001')
iLine2, = plt.loglog(r,(isoMass(r, 0.001)/M_e), '--', label=r'PIM $\rm \alpha_v$ = 0.0001')
#gLine1, = plt.loglog(growthTrack(0.01, 0.01, 0.01, 0.01, 25)[0],growthTrack(0.01, 0.01, 0.01, 0.01, 25)[1], 'b', label=r'($\rm \alpha$, St, $\rm \xi$)=0.01,0.01,0.01')              
gLine1, = plt.loglog(growthTrack(0.01, 0.01, 0.02, 0.01, 25)[0],growthTrack(0.01, 0.01, 0.02, 0.01, 25)[1], 'g', label=r'($\rm \alpha$, St, $\rm \xi$)=0.01,0.01,0.02')
#gLine3, = plt.loglog(growthTrack(0.001, 0.001, 0.01, 0.01, 25)[0],growthTrack(0.001, 0.001, 0.01, 0.01, 25)[1], 'r', label=r'($\rm \alpha$, St, $\rm \xi$)=0.001,0.001,0.01')
#gLine4, = plt.loglog(growthTrack(0.001, 0.001, 0.02, 0.01, 25)[0],growthTrack(0.001, 0.001, 0.02, 0.01, 25)[1], 'y', label=r'($\rm \alpha$, St, $\rm \xi$)=0.001,0.001,0.02')
plt.loglog(growthTrack(0.01, 0.01, 0.02, 0.01, 25)[4],growthTrack(0.01, 0.01, 0.02, 0.01, 25)[5], 'ro')
#plt.loglog(growthTrack(0.01, 0.01, 0.01, 0.01, 5)[0],growthTrack(0.01, 0.01, 0.01, 0.01, 5)[1], 'b')
#plt.loglog(growthTrack(0.01, 0.01, 0.02, 0.01, 5)[0],growthTrack(0.01, 0.01, 0.02, 0.01, 5)[1], 'g')
#plt.loglog(growthTrack(0.001, 0.001, 0.01, 0.01, 5)[0],growthTrack(0.001, 0.001, 0.01, 0.01, 5)[1], 'r')
#plt.loglog(growthTrack(0.001, 0.001, 0.02, 0.01, 5)[0],growthTrack(0.001, 0.001, 0.02, 0.01, 5)[1], 'y')
#plt.loglog(tR_vals,tM_vals, 'ro')

plt.ylabel(r'$M[M_{\rm E}]$')
plt.xlabel(r'$r$[AU]')
plt.legend(fontsize='small', frameon=False, labelspacing=0.2)

#firstLegend = plt.legend(handles=[gLine1, gLine2, gLine3, gLine4], loc=3, fontsize='small', frameon=False, labelspacing=0.01)
firstLegend = plt.legend(handles=[gLine1], loc=3, fontsize='small', frameon=False, labelspacing=0.2)
ax = plt.gca().add_artist(firstLegend)
plt.legend(handles=[iLine1, iLine2], fontsize='small', loc=1, frameon=False, labelspacing=0.1)
plt.plot([0,100], [10,10], ':', color = 'k')
plt.plot([0,100], [20,20], ':', color = 'k')
plt.plot([0,100], [318,318], ':', color = 'k')
#plt.plot([5.2,5.2], [0,3000], ':', color = 'k')
plt.savefig('New6_Step_Int_growthTracks_fullmig.png', dpi=500)
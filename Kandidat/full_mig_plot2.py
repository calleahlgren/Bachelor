#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 00:23:35 2018

@author: calle
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

from full_mig_test import *

r = np.linspace(100, 0.1, 25000)

plt.ylim(0.01, 1000)
plt.xlim(0.1, 100)

plt.loglog(r,(isoMass(r, 0.01)/M_e))
plt.loglog(r,(isoMass(r, 0.001)/M_e))
growthTrack(0.01, 0.01, 0.01, 0.01, 10**-8, 25)
growthTrack(0.01, 0.01, 0.02, 0.01, 10**-8, 25)
growthTrack(0.001, 0.001, 0.01, 0.01, 10**-8, 25)
growthTrack(0.001, 0.001, 0.02, 0.01, 10**-8, 25)
growthTrack(0.01, 0.01, 0.01, 0.01, 10**-8, 5)
growthTrack(0.01, 0.01, 0.02, 0.01, 10**-8, 5)
growthTrack(0.001, 0.001, 0.01, 0.01, 10**-8, 5)
growthTrack(0.001, 0.001, 0.02, 0.01, 10**-8, 5)

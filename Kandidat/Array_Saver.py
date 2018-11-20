#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 18:36:32 2018

@author: calle
"""

#from tempfile import TemporaryFile
import numpy as np
#from Mig_Functions import *
#from Growth_Maps import *

#outfile = TemporaryFile()

#x = np.arange(10)
#np.save("Radial_Pos_Data", Z2)

x = np.load("Mass_Data.npy")
y = np.load("Radial_Pos_Data.npy")

#outfile.seek(0)
#np.load(outfile)
#np.save(outfile, Z1)
#print(Z1)
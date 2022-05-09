#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:55:35 2022

@author: wljw75
"""

"""
Iterates through varying step sizes in r, s.t we can test we have sufficient 
radial bins!
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/wljw75/Documents/phd/KYAGNSED/src')

from kyagnsed import kyagnsed


#Params for testing
M = 1e8
D = 200
log_mdot = -1.2
a = 0.7
cos_inc = 0.5
kte_h =100
kte_w = 0.2
gamma_h = 1.7
gamma_w = 2.4
r_h = 25
r_w = 100
l_rout = -1
fcol = 1
hmax = 10
rep = 0
z = 0



dr_dexs = np.array([120, 110, 100, 90, 80, 70, 60]) #dx to test with
for i in range(len(dr_dexs)):
    agn = kyagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h,
                   gamma_w, r_h, r_w, l_rout, fcol, hmax, rep, z)
    
    agn._change_rBins(dr_dexs[i])
    Lnu_dex = agn.totSpec_std()
    
    if i == 0:
        Lnu_max = Lnu_dex
    
    plt.plot(agn.nu_grid, 100 * (Lnu_max - Lnu_dex)/Lnu_max, label=str(dr_dexs[i]) + ' dex')

plt.xlabel('frequency,   Hz')
plt.ylabel('% Difference')
plt.xscale('log')
plt.legend(frameon=False)
plt.show()
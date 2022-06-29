#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:20:15 2022

@author: wljw75
"""

"""
Varies spacing on energy grid to test convergence
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/wljw75/Documents/phd/KYAGNSED/src')

from relagn import relagnsed


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


def interp_spec(Ei, Emod, flxs):
    """
    Linearly interpolates spectra between energy grids
    Needed for comparing results!

    Parameters
    ----------
    Ei : float
        Point on new energy grid.
    Emod : 1D-array
        Old energy grid
    flxs : 1D-array
        Model spectrum

    Returns
    -------
    fi : float
        Flux at new grid point.

    """
    idx_1 = np.abs(Ei - Emod).argmin()
    
    
    if Ei - Emod[idx_1] > 0:
        if Emod[idx_1] != Emod[-1]: #ensuring we dont fall off array
            E1 = Emod[idx_1]
            E2 = Emod[idx_1 + 1]
            f1 = flxs[idx_1]
            f2 = flxs[idx_1 + 1]
        
        else:
            E1 = Emod[idx_1 - 1]
            E2 = Emod[idx_1]
            f1 = flxs[idx_1 -1]
            f2 = flxs[idx_1]
        
        df_dE = (f2 - f1)/(E2 - E1)
        fi = df_dE * (Ei - E1) + f1
    
    elif Ei - Emod[idx_1] < 0:
        if Emod[idx_1] != Emod[0]:
            E1 = Emod[idx_1 - 1]
            E2 = Emod[idx_1]
            f1 = flxs[idx_1 -1]
            f2 = flxs[idx_1]
        
        else:
            E1 = Emod[idx_1]
            E2 = Emod[idx_1 + 1]
            f1 = flxs[idx_1]
            f2 = flxs[idx_1 + 1]
            
        df_dE = (f2 - f1)/(E2 - E1)
        fi = df_dE * (Ei - E1) + f1
    
    else:
        fi = flxs[idx_1]
        
    return fi


nEs = np.array([5000, 4000, 3000, 2000, 1000, 500])
for i in range(len(nEs)):
    agn = relagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h,
                   gamma_w, r_h, r_w, l_rout, fcol, hmax, rep, z)
    
    agn.new_ear(np.geomspace(1e-4, 1e4, nEs[i]))
    Lnu_dex = agn.totSpec_std()
    
    if i == 0:
        Lnu_max = Lnu_dex
        Lnu_spec = Lnu_max
        nu_mainGrid = agn.nu_grid
    
    elif i != 0:
        Lnu_spec = np.array([])
        for j in range(len(nu_mainGrid)):
            Lj = interp_spec(nu_mainGrid[j], agn.nu_grid, Lnu_dex)
            Lnu_spec = np.append(Lnu_spec, [Lj])
    
    plt.plot(nu_mainGrid, 100 * (Lnu_max - Lnu_spec)/Lnu_max, label=str(nEs[i]) + ' bins')

plt.xlabel('frequency,   Hz')
plt.ylabel('% Difference')
plt.xscale('log')
plt.ylim(-10, 10)
plt.legend(frameon=False)
plt.show()
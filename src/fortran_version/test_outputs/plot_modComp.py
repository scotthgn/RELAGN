#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 17:25:07 2022

@author: wljw75

plots xspec output and compares to python version for same input params
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/wljw75/Documents/phd/RELAGN/src')

from relagn import relagnsed


#params
M = 10
dist = 1e-3
log_mdot = -1
astar = [0, 0.5, 0.998]
cos_inc = 0.8
kTe_hot = 100
kTe_warm = 0.2
gamma_hot = 1.7
gamma_warm = 2.7
r_hot = 10
r_warm = 20
log_rout = 7
fcol = 1
hmax = 10
z = 0


Es0, fs0 = np.loadtxt('modtest_0spin.qdp', skiprows=3, usecols=(0, 2), unpack=True)
Es05, fs05 = np.loadtxt('modtest_05spin.qdp', skiprows=3, usecols=(0, 2), unpack=True)
Esm, fsm = np.loadtxt('modtest_maxspin.qdp', skiprows=3, usecols=(0, 2), unpack=True)


plt.loglog(Es0, fs0, color='red', label='0 spin')
plt.loglog(Es05, fs05, color='green', label='0.5 spin')
plt.loglog(Esm, fsm, color='blue', label='0.998 spin')

for a in astar:
    ragn = relagnsed(M, dist, log_mdot, a, cos_inc, kTe_hot, kTe_warm, gamma_hot,
                     gamma_warm, r_hot, r_warm, log_rout, fcol, hmax, 0, z)
    
    ragn.set_counts()
    ragn.set_flux()
    
    fspy = ragn.totSpec_rel()
    espy = ragn.E_obs
    
    if a == 0.998:
        plt.loglog(espy, espy**2 * fspy, ls='dashed', color='k', label='Python version')
    else:
        plt.loglog(espy, espy**2 * fspy, ls='dashed', color='k')
    

plt.ylim(1, 600)

plt.xlabel('Energy   (keV)')
plt.ylabel(r'EF(E)   keV$^{2}$ (Photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$)')
plt.legend(frameon=False)
plt.show()
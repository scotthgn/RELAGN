#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 17:18:34 2022

@author: wljw75

Plots the xspec output comparing to kerrbb
"""

import numpy as np
import matplotlib.pyplot as plt


Es0, fs0, fsk0 = np.loadtxt('comp2kerrbb_0spin.qdp', skiprows=3, usecols=(0, 3, 4),
                            unpack=True)
Es05, fs05, fsk05 = np.loadtxt('comp2kerrbb_05spin.qdp', skiprows=3, usecols=(0, 3, 4),
                            unpack=True)
Esm, fsm, fskm = np.loadtxt('comp2kerrbb_maxspin.qdp', skiprows=3, usecols=(0, 3, 4),
                            unpack=True)



plt.loglog(Es0, fs0, color='red', label='0 spin')
plt.loglog(Es0, fsk0, ls='dashed', color='k')

plt.loglog(Es05, fs05, color='green', label='0.5 spin')
plt.loglog(Es05, fsk05, ls='dashed', color='k')

plt.loglog(Esm, fsm, color='blue', label='0.998 spin')
plt.loglog(Esm, fskm, ls='dashed', color='k', label='kerrbb')

plt.ylim(1, 600)
plt.xlabel('Energy   (keV)')
plt.ylabel(r'EF(E)   keV$^{2}$ (Photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$)')
plt.legend(frameon=False)
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 15:47:19 2022

@author: wljw75

Simply plots the SED you want easily - so we can compare to fortran version
"""

import matplotlib.pyplot as plt

import sys
sys.path.append('/home/wljw75/Documents/phd/RELAGN/src')

from relagn import relagnsed


M = 10
dist = 1e-3
log_mdot = -1
a = 0.998
cos_inc = 0.5   
kTe_hot = 100
kTe_warm = 0.2
gamma_hot = 1.7
gamma_warm = 2.7
r_hot = 1
r_warm = 2
log_rout = 7
fcol = 1
h_max = 10
reprocess = 0
z = 0

ragn = relagnsed(M, dist, log_mdot, a, cos_inc, kTe_hot, kTe_warm, gamma_hot,
                 gamma_warm, r_hot, r_warm, log_rout, fcol, h_max, reprocess, z)

ragn.set_counts()
ragn.set_flux()
Es = ragn.E_obs

ftot = ragn.totSpec_rel()
fd = ragn.Lnu_disc_rel
fw = ragn.Lnu_warm_rel
fh = ragn.Lnu_hot_rel

ftot_nr = ragn.totSpec_std()
fdnr = ragn.Lnu_disc_norel

plt.loglog(Es, Es**2 * fd)
plt.loglog(Es, Es**2 * fdnr)
plt.xlim(1e-2, 1e2)
plt.ylim(10, 1e3)
plt.show()
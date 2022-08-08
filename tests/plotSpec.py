#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 15:47:19 2022

@author: wljw75

Simply plots the SED you want easily - so we can compare to fortran version
"""

import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

import sys
sys.path.append('/home/wljw75/Documents/phd/RELAGN/src')

from relagn import relagnsed


M = 1e5
dist = 1
log_mdot = -1
a = 0
cos_inc = 0.5   
kTe_hot = 100
kTe_warm = 0.2
gamma_hot = 1.7
gamma_warm = 2.7
r_hot = 10
r_warm = 20
log_rout = -1
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
fwnr = ragn.Lnu_warm_norel

#plt.loglog(Es, Es**2 * fd)
#plt.loglog(Es, Es**2 * fdnr)
plt.loglog(Es, Es**2 * (fd + fw))
plt.xlim(1e-4, 10)
plt.ylim(1e-2, 3)
plt.show()

sigma_sb = 5.670367e-5 #erg/s/cm^2/K^4
def normw_ann(r, dr):
    t4 = ragn.calc_Tnt(r)
    d = (dist*u.Mpc).to(u.cm).value
    
    norm = sigma_sb * t4
    norm = norm * 4*np.pi*r*dr * (ragn.Rg*100)**2
    norm = norm * cos_inc/0.5
    norm = norm/(4*np.pi * d**2)
    return norm

print(normw_ann(19.318726578496907, 1.3393401692638527))

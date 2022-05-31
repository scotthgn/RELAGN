#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 10:29:44 2022

@author: wljw75
"""

"""
Creates spectrum with relconv - and slow on without to check is sensible...
"""

import numpy as np
import matplotlib.pyplot as plt
import xspec
import astropy.units as u

import sys
sys.path.append('/home/wljw75/Documents/phd/KYAGNSED/src/')

from kyagnsed import kyagnsed
import time

ts = time.time()


#Params for testing
M = 10
D = 1e-2
log_mdot = -1.2
a = 0.7
cos_inc = 0.9
kte_h =100
kte_w = 0.2
gamma_h = 1.7
gamma_w = 2.4
r_h = -1
r_w = -1
l_rout = -1
fcol = 1
hmax = 10
rep = 0
z = 0


fig = plt.figure(figsize=(6, 5))
ax1 = fig.add_subplot(111)

#Calculating my model components
myagn = kyagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h, gamma_w,
                  r_h, r_w, l_rout, fcol, hmax, rep, z)

myagn.set_counts() #sets cgs units
myagn.set_flux() #units of flux
nu = myagn.E_obs

#extracting rel components
ftot_r = myagn.totSpec_rel()
ax1.loglog(nu, nu*ftot_r, color='red', label='relativistic')

try:
    fd_r = myagn.Lnu_disc_rel
    ax1.loglog(nu, nu*fd_r, color='red')
except:
    pass

try:
    fw_r = myagn.Lnu_warm_rel
    ax1.loglog(nu, nu*fw_r, color='green')
except:
    pass

try:
    fh_r = myagn.Lnu_hot_rel
    ax1.loglog(nu, nu*fh_r, color='blue')
except:
    pass


#Extracting non-rel components
ftot_n = myagn.totSpec_std()
ax1.loglog(nu, nu*ftot_n, color='k', ls='-.', label='non-relativistic')

try:
    fd_n = myagn.Lnu_disc_norel
    ax1.loglog(nu, nu*fd_n, color='k', ls='-.')
except:
    pass

try:
    fw_n = myagn.Lnu_warm_norel
    ax1.loglog(nu, nu*fw_n, color='green', ls='-.')
except:
    pass

try:
    fh_n = myagn.Lnu_hot_norel
    ax1.loglog(nu, nu*fh_n, color='blue', ls='-.')
except:
    pass

ax1.legend(frameon=False)
#ax1.set_ylim(1e-8, 1e-5)
ax1.set_ylim(1, 1e1)
ax1.set_xlim(1e-3, 10)
#ax1.set_xlim(3e13, 1e21)
ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel(r'EF(E)   keV$^{2}$ (Photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$)')

"""
def to_energy(nu):
    return (nu*u.Hz).to(u.keV, equivalencies=u.spectral()).value

ax2 = ax1.secondary_xaxis('top', functions=(to_energy, to_energy))
ax2.set_xlabel('Energy (keV)')
"""

tf = time.time()
print('Runtime = {}s'.format(tf - ts))

plt.tight_layout()
plt.show()


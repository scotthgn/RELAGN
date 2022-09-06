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
sys.path.append('/home/wljw75/Documents/phd/RELAGN/src/')

from relagn import relagnsed
import time

ts = time.time()


#Params for testing
M = 2e8
D = 200
log_mdot = -1.16700
a = 0.773263
cos_inc = 0.9
kte_h =100
kte_w = 0.286181
gamma_h = 1.95805
gamma_w = 2.69997
r_h = 17.6806
r_w = 1000
l_rout = -1
fcol = 1
hmax = 10
rep = 1
z = 0.045

fig = plt.figure(figsize=(6, 5))
ax1 = fig.add_subplot(111)

#Calculating my model components
myagn = relagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h, gamma_w,
                  r_h, r_w, l_rout, fcol, hmax, rep, z)

myagn.set_counts() #sets cgs units
myagn.set_flux() #units of flux
nu = myagn.E_obs

#extracting rel components
ftot_r = myagn.totSpec_rel()
ax1.loglog(nu, nu**2 * ftot_r, color='red', label='relativistic')

try:
    fd_r = myagn.Lnu_disc_rel
    ax1.loglog(nu, nu**2 * fd_r, color='red')
except:
    pass

try:
    fw_r = myagn.Lnu_warm_rel
    ax1.loglog(nu, nu**2 * fw_r, color='green')
except:
    pass

print(myagn.Lnu_hot_rel)
try:
    fh_r = myagn.Lnu_hot_rel
    ax1.loglog(nu, nu**2 * fh_r, color='blue')
except:
    pass


#Extracting non-rel components
ftot_n = myagn.totSpec_std()
ax1.loglog(nu, nu**2 * ftot_n, color='k', ls='-.', label='non-relativistic')

try:
    fd_n = myagn.Lnu_disc_norel
    ax1.loglog(nu, nu**2 * fd_n, color='k', ls='-.')
except:
    pass

try:
    fw_n = myagn.Lnu_warm_norel
    ax1.loglog(nu, nu**2 * fw_n, color='green', ls='-.')
except:
    pass

try:
    fh_n = myagn.Lnu_hot_norel
    ax1.loglog(nu, nu**2 * fh_n, color='blue', ls='-.')
except:
    pass

ax1.legend(frameon=False)
#ax1.set_ylim(1e-10, 1e-8)
ax1.set_ylim(1e-3, 1e-1)
ax1.set_xlim(1e-3, 1e3)
#ax1.set_xlim(3e13, 1e20)
ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel(r'EF(E)   keV$^{2}$ (Photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$)')


def to_energy(nu):
    return (nu*u.Hz).to(u.keV, equivalencies=u.spectral()).value

ax2 = ax1.secondary_xaxis('top', functions=(to_energy, to_energy))
ax2.set_xlabel('Energy (keV)')


tf = time.time()
print('Runtime = {}s'.format(tf - ts))

plt.tight_layout()
plt.show()


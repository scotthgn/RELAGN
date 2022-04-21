#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:48:44 2022

@author: wljw75
"""

"""
Tests relativistic disc agianst kerrbb. For disc extending from r_isco to
infinity - these two models should be identical
"""

import numpy as np
import matplotlib.pyplot as plt
import xspec
import astropy.units as u

import sys
sys.path.append('/home/wljw75/Documents/phd/RELAGNSED/src/')

from relagnsed import relagnsed


#Params for testing
M = 10
D = 1
log_mdot = -1
a = 0
cos_inc = 0.5
kte_h = 100
kte_w = 1
gamma_h = 1.8
gamma_w = 2.2
r_h = -1   
r_w = -1
l_rout = 5
fcol = 1
hmax = 10
rep = 0
z = 0

#Calculating my model components
myagn = relagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h, gamma_w,
                  r_h, r_w, l_rout, fcol, hmax, rep, z)

myagn.set_cgs()
myagn.set_flux()

nus = myagn.nu_obs
fs = myagn.totSpec_rel()
fs_n = myagn.totSpec_std() #for comparison!

#Now douing kerrbb
xspec.Xset.chatter = 0 
xspec.AllData.dummyrsp(myagn.Emin, myagn.Emax, myagn.numE)

Mdot = myagn.mdot * myagn.Mdot_edd * 1e3
Dkpc = (myagn.D * u.Mpc).to(u.kpc).value
kerpars = (0, a, np.rad2deg(myagn.inc), M, Mdot/1e18, Dkpc, 1, rep, 0)
xspec.Model('kerrbb', setPars=kerpars)

xspec.Plot('model')
kes = np.array(xspec.Plot.x())
kph = np.array(xspec.Plot.model())
 
xspec.AllModels.clear()
    
fs_kerr = kph * kes #keV Photons/s/cm^2/keV
fs_kerr = (fs_kerr * u.keV/u.s/u.keV).to(u.erg/u.s/u.Hz, equivalencies=u.spectral())
nu_kerr = (kes * u.keV).to(u.Hz, equivalencies=u.spectral())


#Plotting
fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(111)

ax1.loglog(nus, nus*fs, color='red', label='relconv*agnsed_disc')
ax1.loglog(nu_kerr, nu_kerr*fs_kerr, color='k', ls='-.', label='kerrbb')

ax1.loglog(nus, nus*fs_n, color='blue', label='agnsed_disc')

ax1.legend(frameon=False)

ax1.set_ylim(1e-14, 1e-12)
ax1.set_xlim(3e15, 1e19)

ax1.set_xlabel(r'Frequency, $\nu$   (Hz)')
ax1.set_ylabel(r'$\nu F_{\nu}$   (ergs cm$^{-2}$ s$^{-1}$)')


def to_energy(nu):
    return (nu*u.Hz).to(u.keV, equivalencies=u.spectral()).value

ax2 = ax1.secondary_xaxis('top', functions=(to_energy, to_energy))
ax2.set_xlabel('Energy (keV)')

plt.tight_layout()
plt.show()
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:18:49 2022

@author: wljw75
"""

"""
Tests the intrinsic non-relativsic calculations agains AGNSED, as these
should be identical
"""

import numpy as np
import matplotlib.pyplot as plt
import xspec
import astropy.units as u

import sys
sys.path.append('/home/wljw75/Documents/phd/RELAGN/src/')

from relagn import relagnsed


#Params for testing
M = 2e8
D = 200
log_mdot = -1.27096
a = 0
cos_inc = 0.9
kte_h =100
kte_w = 0.248287
gamma_h = 1.97568
gamma_w = 2.65004
r_h = 24.4264
r_w = 1e3
l_rout = -1
fcol = 1
hmax = 10
rep = 0
z = 0


#Calculating my model components
myagn = relagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h, gamma_w,
                  r_h, r_w, l_rout, fcol, hmax, rep, z)

myagn.set_counts() #sets cgs units
myagn.set_flux()

nus = myagn.nu_obs
Es = myagn.E_obs
ftot = myagn.totSpec_std() 

fd = myagn.Lnu_disc_norel
fw = myagn.Lnu_warm_norel
fh = myagn.Lnu_hot_norel


#Now calculating AGNSED version
xspec.Xset.chatter = 0 
xspec.AllData.dummyrsp(myagn.Emin, myagn.Emax, myagn.numE)

agnpars = (M, D, log_mdot, a, cos_inc, kte_h, 0.2, gamma_h, gamma_w, r_h,
           500, l_rout, hmax, rep, z)

mx1 = xspec.Model('agnsed', setPars=agnpars)
mx1.setPars({7:str(kte_w) + ' 0.01 0.1 0.1 2 2'})
mx1.setPars({11:str(r_w) + ' 1 6 6 1e3 1e3'})

xspec.Plot('model')
es_agn = np.array(xspec.Plot.x())
ph = np.array(xspec.Plot.model())

#fs_agn = ph * es
#fs_agn = (fs_agn * u.keV/u.s/u.keV).to(u.erg/u.s/u.Hz, equivalencies=u.spectral()).value

#nu_agn = (es * u.keV).to(u.Hz, equivalencies=u.spectral())

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(111)

ax1.loglog(Es, Es**2 * fd, color='red', ls='-.', label='Disc')
ax1.loglog(Es, Es**2 * fw, color='green', ls='-.', label='Warm Comp.')
ax1.loglog(Es, Es**2 * fh, color='blue', ls='-.', label='Hot Comp.')
ax1.loglog(Es, Es**2 * ftot, color='k', label='Total')

ax1.loglog(es_agn, es_agn**2 * ph, color='gray', ls='-.', label='XSPEC (agnsed)')

ax1.legend(frameon=False)

ax1.set_ylim(1e-3, 0.1)
#ax1.set_ylim(1e-10, 1e-7)
ax1.set_xlim(1e-4, 1e4)

ax1.set_xlabel(r'Frequency, $\nu$   (Hz)')
ax1.set_ylabel(r'$\nu F_{\nu}$   (ergs cm$^{-2}$ s$^{-1}$)')


def to_energy(nu):
    return (nu*u.Hz).to(u.keV, equivalencies=u.spectral()).value

ax2 = ax1.secondary_xaxis('top', functions=(to_energy, to_energy))
ax2.set_xlabel('Energy (keV)')

ax1.axvline(4e17, ls='dotted', color='green')

plt.tight_layout()
plt.show()


plt.plot(nus, 100 * (ftot-fs_agn)/ftot)
plt.xlim(1e13, 1e19)
plt.ylim(-2.5, 2.5)
plt.xscale('log')
plt.axvline(4e17, ls='dotted', color='green')
plt.axhline(0, ls='-.', color='k')
plt.ylabel('% Difference')
plt.show()  

plt.loglog(nus, abs(nus/nu_agn.value))
plt.xlim(1e13, 1e19)
plt.show()

print(np.trapz(ftot, nus), np.trapz(fs_agn, nu_agn))
print(myagn.r_sg)

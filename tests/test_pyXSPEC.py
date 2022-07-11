#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 13:08:39 2022

@author: wljw75
"""

"""
Tests that the import into pyXSPEC has worked as it should
"""

import numpy as np
import matplotlib.pyplot as plt
import xspec

import sys
sys.path.append('/home/wljw75/Documents/phd/RELAGN/src/')

from relagn import relagnsed
import load_2xspec


M = 10
D = 1
log_mdot = -1
a = 0
cos_inc = 0.5
kte_h = 100
kte_w = 1
gamma_h = 1.8
gamma_w = 2.2
r_h = 10
r_w = 20
l_rout = 7
fcol = 1
hmax = 10
rep = 0
z = 0

ragn = relagnsed(M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h, gamma_w,
                 r_h, r_w, l_rout, fcol, hmax, rep, z)

ragn.set_counts()
ragn.set_flux()

Eagn = ragn.Egrid
ph_agn = ragn.totSpec_rel()


#Now for pyXSPEC verion
pars = (M, D, log_mdot, a, cos_inc, kte_h, kte_w, gamma_h, gamma_w, r_h, r_w,
        l_rout, fcol, hmax, rep)
xspec.AllData.dummyrsp(1, 10, 500)
xspec.Model('relagn', setPars=pars)

xspec.Plot.device = '/null'
xspec.Plot('model')
es = np.array(xspec.Plot.x())
fs = np.array(xspec.Plot.model())

plt.loglog(Eagn, Eagn**2 * ph_agn)
plt.loglog(es, es**2 * fs)
plt.ylim(1e-6, 1e-3)
plt.show()

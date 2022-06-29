#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 12:04:18 2022

@author: wljw75
"""
"""
Calls relagnsed and saves output to a .dat file
Yes, this is a faf, but unfortuanetly necessary to avoid segmentation fault
in xspec...
"""

import numpy as np
from relagn import relagnsed
import os
import matplotlib.pyplot as plt


wdir = os.getcwd()

params = np.loadtxt(wdir + '/relpars.dat')
M = params[0]
dist = params[1]
log_mdot = params[2]
a = params[3]
cos_inc = params[4]
kTe_hot = params[5]
kTe_warm = params[6]
gamma_hot = params[7]
gamma_warm = params[8]
r_hot = params[9]
r_warm = params[10]
log_rout = params[11]
fcol = params[12]
h_max = params[13]
reprocess = params[14]
z = params[15] #fixing redshift at 0 for now - change later

ragn = relagnsed(M, dist, log_mdot, a, cos_inc, abs(kTe_hot), abs(kTe_warm), 
                 gamma_hot, abs(gamma_warm), r_hot, r_warm, log_rout, fcol, h_max, reprocess, z)

ear = np.loadtxt(wdir + '/ear.dat')
Els = np.array(ear[:-1])
Ers = np.array(ear[1:])
dEs = Ers - Els
Emid = Els + 0.5*dEs #evaluate model at center of bin


#Checking wheather to extend default energy grid
if min(ear) < 1e-4:
    new_emin = min(ear)
else:
    new_emin = 1e-4

if max(ear) > 1e4:
    new_emax = max(ear)
else:
    new_emax = 1e4


ragn.new_ear(np.geomspace(new_emin, new_emax, 1000)) #Ensuring sufficient energy coverage

ragn.set_counts()
ragn.set_flux() #getting output as photons/s/cm^2/keV

Emod = ragn.E_obs

flxs_all = ragn.totSpec_rel() #total spectrum

#Checking switching params
if kTe_hot < 0:
    flxs = ragn.Lnu_hot_rel
elif kTe_warm < 0:
    flxs = ragn.Lnu_warm_rel
elif gamma_warm < 0:     
    flxs = ragn.Lnu_disc_rel
else:
    flxs = flxs_all

idx_neg = np.argwhere(flxs < 0)
flxs[idx_neg] = 0.


#Extracting part of the spectrum of interest
idx_low = np.argwhere(Emod < min(ear))
Emod = np.delete(Emod, idx_low)
flxs = np.delete(flxs, idx_low)

idx_high = np.argwhere(Emod > max(ear))
Emod = np.delete(Emod, idx_high)
flxs = np.delete(flxs, idx_high)


#re-casting onto original grid through linear interpolation
def interp_spec(Ei):
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

flx_spec = np.array([])
for i in range(len(Emid)):
    fnew = interp_spec(Emid[i])
    flx_spec = np.append(flx_spec, [dEs[i] * fnew])
    
    
if any(np.isnan(flx_spec)):
    print(params)

np.savetxt('relagn_fluxs.dat', flx_spec)
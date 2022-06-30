#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:58:37 2022

@author: wljw75

Contains the functions necessary for loading relagn to xspec. Includes both 
relativistic and non-relativistc functions
"""

import numpy as np
import xspec
from relagn import relagnsed



def interp_spec(Ei, Emod, flxs):
    """
    Linearly interpolates spectra between energy grids
    Needed for tunring into xspec model...

    Parameters
    ----------
    Ei : float
        Point on new energy grid.
    Emod : 1D-array
        Old energy grid
    flxs : 1D-array
        Model spectrum

    Returns
    -------
    fi : float
        Flux at new grid point.

    """
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



def relagn(ear, params, flx):
    
    #Read params
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
    z = params[15] 
    
    #Initiate model
    rsed = relagnsed(M, dist, log_mdot, a, cos_inc, kTe_hot, kTe_warm, 
                 gamma_hot, gamma_warm, r_hot, r_warm, log_rout, fcol, 
                 h_max, reprocess, z)
        
    #Extracting energy bins
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
        
    rsed.new_ear(np.geomspace(new_emin, new_emax, 1000)) #Ensuring sufficient energy coverage
    
    rsed.set_counts()
    rsed.set_flux() #getting output as photons/s/cm^2/keV
    
    #model
    Emod = rsed.E_obs
    flxs_all = rsed.totSpec_rel() #total spectrum
    
    #Checking switching params
    if kTe_hot < 0:
        flxs = rsed.Lnu_hot_rel
    elif kTe_warm < 0:
        flxs = rsed.Lnu_warm_rel
    elif gamma_warm < 0:     
        flxs = rsed.Lnu_disc_rel
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
    for i in range(len(ear) - 1):
        fnew = interp_spec(Emid[i], Emod, flxs)
        flx[i] = fnew * dEs[i]
        
    
    
#Now making local model
parinf = ('M Msol 10 1 1 1e10 1e10 -1',
          'dist Mpc 1 1e-3 1e-3 1e3 1e3 -1',
          'log_mdot "" -1 -10 -10 2 2 0.01',
          'astar "" 0 0 0 0.998 0.998 0.001',
          'cosi "" 0.5 0.09 0.09 0.998 0.998 0.001',
          'kTe_hot keV 100 10 10 300 300 -1',
          'kTe_warm keV 0.2 0.01 0.01 1 1 0.01',
          'gamma_hot "" 1.7 1.3 1.3 3 3 0.01',
          'gamma_warm "" 2.7 2 2 5 10 0.01',
          'r_hot Rg 10 6 6 500 500 0.01',
          'r_warm Rg 20 6 6 500 500 0.01',
          'logrout self_g -1 -1 -1 7 7 -1',
          'fcol "" 1 1 1 5 5 0.01',
          'hmax Rg 10 6 6 10 10 -1',
          '$reprocess 0',
          'z "" 0 0 0 1 1 -1')
    
xspec.AllModels.addPyMod(relagn, parinf, 'add')

print('relagn loaded succesfully')
print('If you use this in your work please reference Hagen & Done (in prep.)')
print()
print('Call as relagn')
    
    
    
    
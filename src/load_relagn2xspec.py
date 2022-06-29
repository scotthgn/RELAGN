#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 09:29:24 2022

@author: wljw75
"""

"""
Calls kyagnsed for given set of input parameters + energy array and returns
the spectrum. This is to make it considerably easier to load as a local model
in pyXSPEC. Includes the function load_kyagn2xspec() - which calls imports relagnsed
as a local model in xspec
"""

import xspec
import numpy as np
import os

wdir = os.getcwd()
fpath = os.path.dirname(os.path.realpath(__file__))
def relagn(ear, params, flx):    
    np.savetxt('ear.dat', ear)
    np.savetxt('relpars.dat', params)
    os.system('python ' + fpath + '/call_relagn.py')
    
    flxs = np.loadtxt(wdir + '/relagn_fluxs.dat')
    
    os.system('rm ear.dat')
    os.system('rm relpars.dat')
    os.system('rm relagn_fluxs.dat')
    
    for i in range(len(ear) - 1):
        flx[i] = flxs[i] #photons/s/cm^s
        
    
#parameter info for xspec
parinfo_ragn = ('M Msol 10 1 1 1e10 1e10 -1',
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

 
xspec.AllModels.addPyMod(relagn, parinfo_ragn, 'add')
print('relagnsed model loaded succefully')
print('in pyXSPEC call as relagn')
    


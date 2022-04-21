#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 11:46:24 2022

@author: wljw75
"""

"""
Tests the relconv convolution on a simple gaussian with a r^-3 emissivity
law - as this should (hoefully!) give the same as kdblur
"""

import numpy as np
import xspec
import matplotlib.pyplot as plt



xspec.Xset.chatter = 0
xspec.AllModels.lmod('relxill', dirPath='/home/wljw75/relxill/')

"""Defining model parameters"""
r_in = 1.236
r_out = 100
inc = 30
sigma = 0.01
El = 1




"""Testing with kdblur * gauss"""
gauss_pars = (3, r_in, r_out, inc, El, sigma, 1)
gauss_mo = xspec.Model('kdblur*gauss', setPars=gauss_pars)
xspec.AllData.dummyrsp(1e-2, 1e4)

xspec.Plot.device = '/null'
xspec.Plot('model')

kd_es = np.array(xspec.Plot.x())
kd_fs = np.array(xspec.Plot.model()) * 0.01

xspec.AllModels.clear()


"""Testing relconv on xspec gauss"""
#relparameters are: idx1, idx2, r_br, a, inc, r_in, r_out, limb
relpars = (3, 3, 15, 0.998, inc, -1*r_in/1.236, r_out, 1)
all_relgpars = relpars + (El, sigma, 1)

mo_relgauss = xspec.Model('relconv*gauss', setPars=all_relgpars)
xspec.Plot('model')

rx_es = np.array(xspec.Plot.x())
rx_fs = np.array(xspec.Plot.model())






"""
Now testing with my own model and relconv
"""
def emiss(rls, rhs):
    dr = rhs - rls
    rm = rls + 0.5*dr
    eps = (rm/rls[0])**(-3)
    int_area = np.sum(np.pi * (rls + rhs) * (rhs - rls) * eps)
    
    eps *= (2*np.pi*rm*dr)/int_area
    return eps

def relconv_gauss(rin, rout):
    rs = np.geomspace(rin, rout, 1000)
    norm_rs = emiss(rs[:-1], rs[1:])
    #print(norm_rs)
    
    for i in range(len(rs) - 1):
        dr = rs[i + 1] - rs[i]
        rmid = rs[i] + dr/2
        
        #Creating pyxspec funtion to calculate gaussian at annulus
        def my_gauss(es, params, flx):
            El, sigma, norm = params
            norm_r = 2*np.pi *  rmid* dr * (rmid/rin)**(-3)
            
            Els = np.array(es[:-1])
            Ers = np.array(es[1:])

            dEs = Ers - Els
            Emids = Els + dEs/2

            fluxs =  2e-3 * norm_rs[i] * (1/(sigma * np.sqrt(2*np.pi))) * np.exp(
                    (-(Emids - El)**2)/(2 * sigma**2))
            
            for j in range(len(es) - 1):                
                flx[j] = fluxs[j]
        
        
        parinfo = ('El "keV" 1 0.1 0.1 0.1 10 10',
                   'sigma "keV" 0.1 0.01 0.01 0.01 5 5')
        
        xspec.AllModels.addPyMod(my_gauss, parinfo, 'add')  
        
        #convolving annulus gaussian with relconv
        relparam = (3, 3, rs[i]+dr/2, 0.998, inc, -1*rs[i]/1.236, rs[i]+dr, 1)
        #relparam = (10, 0.998, inc, -1*rs[i]/1.236, rs[i]+dr, 1, 2)
        gparam = (El, sigma)
        all_param = relparam + gparam
        
        if i == 0 or i == len(rs)-1:
            print(relparam)
        
        mo = xspec.Model('relconv*my_gauss', setPars=all_param)     
        
        xspec.Plot('model')
        
        g_es = np.array(xspec.Plot.x())
        g_cs = np.array(xspec.Plot.model())
        
        xspec.AllModels.clear()
        
        
        if i == 0:
            all_gs = g_cs
        else:
            all_gs = np.column_stack((all_gs, g_cs))
        
    #returning the total line profile from the entire disc
    gs_tot = np.sum(all_gs, axis=-1)
    return gs_tot, g_es

fs_g, es_g = relconv_gauss(r_in, r_out)

ftot = np.trapz(fs_g, es_g)
print(ftot)
print(ftot * 2*np.pi*1.236**2)

"""Testing with relconv on single gaussian"""

def m_gauss(es, param, flx):
    El, sigma, norm = param

    
    Els = np.array(es[:-1])
    Ers = np.array(es[1:])
    
    dEs = Ers - Els
    Emids = Els + dEs/2

    fluxs = 0.002 * (1/(sigma * np.sqrt(2*np.pi))) * np.exp(
        (-(Emids - El)**2)/(2 * sigma**2))
            
    for j in range(len(es) - 1):                
        flx[j] = fluxs[j]


parinfo = ('El "keV" 1 0.1 0.1 0.1 10 10',
           'sigma "keV" 0.1 0.01 0.01 0.01 5 5')
        
xspec.AllModels.addPyMod(m_gauss, parinfo, 'add')  
        

relmo = xspec.Model('relconv * m_gauss', setPars=all_relgpars)

xspec.Plot('model')

es_rel = np.array(xspec.Plot.x())
fs_rel = np.array(xspec.Plot.model())








plt.loglog(kd_es, kd_es**2 * kd_fs, label='kdblur*gauss')
plt.loglog(rx_es, rx_es**2 * rx_fs, label='relconv*gauss')
plt.loglog(es_g, es_g**2 * fs_g, label='relconv*my_gauss')
#plt.loglog(es_rel, es_rel**2 * fs_rel, label='relconv*single_my_gauss')

plt.ylim(1e-3, 1e1)
plt.xlim(0.1, 1.5)
plt.xlabel('Energy (keV)')
plt.ylabel('E F(E)')
plt.legend(frameon=False)
plt.title('Physical norm off')

plt.show()

    

    
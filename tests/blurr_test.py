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

"""Defining model parameters"""
r_in = 1.236
r_out = 100
inc = 30
sigma = 0.01
El = 1




"""Testing with kdblur * gauss"""
gauss_pars = (3, r_in, r_out, inc, El, sigma, 1)
gauss_mo = xspec.Model('kdblur*gauss', setPars=gauss_pars)
#gauss_mo = xspec.Model('gauss', setPars=(El, sigma, 1))
xspec.AllData.dummyrsp(0.1, 2, 500)

xspec.Plot.device = '/null'
xspec.Plot('model')

kd_es = np.array(xspec.Plot.x())
kd_fs = np.array(xspec.Plot.model()) 

xspec.AllModels.clear()


"""Testing kyconv on xspec gauss"""
#relparameters are: a, inc, r_in, ms (set 0 to integrate from rin), rout, po_idx1, po_idx2, rb, z, nE, norm
#nE needs to be set to same or more as response matrix
#set norm to -1 to have proper physical normalisation...
relpars = (0.998, inc, r_in, 0, r_out, 3, 3, 10, 0, 0, 500, -1) 
all_relgpars = relpars + (El, sigma, 1)

mo_relgauss = xspec.Model('kyconv*gauss', setPars=all_relgpars)
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
    rs = np.geomspace(rin, rout, 200)
    norm_rs = emiss(rs[:-1], rs[1:])
    #print(norm_rs)
    
    for i in range(len(rs) - 1):
        dr = rs[i + 1] - rs[i]
        rmid = rs[i] + dr/2
        
        #Creating pyxspec funtion to calculate gaussian at annulus
        def my_gauss(es, params, flx):
            El, sigma, norm = params
            norm_r = rmid**(-3)
            #norm_r = ((rin**(-2) - rout**(-2))/(rs[i]**(-2) - rs[i+1]**(-2))) * rmid**(-3)
            
            Els = np.array(es[:-1])
            Ers = np.array(es[1:])

            dEs = Ers - Els
            Emids = Els + dEs/2

            fluxs =  norm_r*(1/(sigma * np.sqrt(2*np.pi))) * np.exp(
                    (-(Emids - El)**2)/(2 * sigma**2))
            #print('flx_mg = {}'.format(np.trapz(fluxs, Emids)))
            for j in range(len(es) - 1):                
                flx[j] = dEs[j]*fluxs[j]
        
        
        parinfo = ('El "keV" 1 0.1 0.1 0.1 10 10',
                   'sigma "keV" 0.1 0.01 0.01 0.01 5 5')
        
        xspec.AllModels.addPyMod(my_gauss, parinfo, 'add')  
        
        #convolving annulus gaussian with relconv
        relparam = (0.998, inc, rs[i], 0, rs[i+1], 0, 0, rmid, 0, 0, 500, -1)
        #relparam = (10, 0.998, inc, -1*rs[i]/1.236, rs[i]+dr, 1, 2)
        gparam = (El, sigma, 1)
        all_param = relparam + gparam
        
        if i == 0 or i == len(rs)-1:
            print(relparam)
        
        mo = xspec.Model('kyconv*my_gauss', setPars=all_param)  
        #mo = xspec.Model('my_gauss', setPars=gparam)
        
        xspec.Plot('model')
        
        g_es = np.array(xspec.Plot.x())
        g_cs = np.array(xspec.Plot.model())
        #print('flx_ky = {}'.format(np.trapz(g_cs, g_es)))
        #print()
        xspec.AllModels.clear()
        
        plt.loglog(g_es, g_es**2 * g_cs, ls='-.')
        
        
        if i == 0:
            all_gs = g_cs
        else:
            all_gs = np.column_stack((all_gs, g_cs))
        
    #returning the total line profile from the entire disc
    gs_tot = np.sum(all_gs, axis=-1)
    plt.loglog(g_es, g_es**2 * gs_tot, color='k')
    plt.ylim(1e-4, 1e1)
    plt.xlim(0.1, 1.5)
    plt.show()
    return gs_tot, g_es

fs_g, es_g = relconv_gauss(r_in, r_out)

ftot = np.trapz(fs_g, es_g)
print(ftot)
print(ftot * 2*np.pi*1.236**2)

"""Testing with relconv on single gaussian"""
"""
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

"""






plt.loglog(kd_es, kd_es**2 * kd_fs, label='kdblur*gauss')
plt.loglog(rx_es, rx_es**2 * rx_fs, label='kyconv*gauss')
plt.loglog(es_g, es_g**2 * fs_g, label='kyconv*my_gauss')
#plt.loglog(es_rel, es_rel**2 * fs_rel, label='relconv*single_my_gauss')

plt.ylim(1e-3, 1e1)
plt.xlim(0.1, 1.5)
plt.xlabel('Energy (keV)')
plt.ylabel('E F(E)')
plt.legend(frameon=False)

plt.show()

    

    
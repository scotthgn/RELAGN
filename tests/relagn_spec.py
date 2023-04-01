#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 15:20:33 2022

@author: Scott Hagen

code to play around with the relagn model - test to see if it gives spectra,
unit handling, etc
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import sys
mdir = os.path.abspath(__file__)
mdir = mdir.replace('/tests/relagn_spec.py', '/src/python_version')
sys.path.append(mdir)

from relagn import relagn


def main():
    #Test pars
    M=1e8
    dist=100
    log_mdot=-1
    a=0.
    cos_inc=0.5
    kTe_hot=100
    kTe_warm=0.2
    gamma_hot=1.7
    gamma_warm=2.7
    r_hot=10
    r_warm=20
    log_rout=-1
    fcol=1
    h_max=10
    z=0
    
    userel = False
    ragn = relagn(M, dist, log_mdot, a, cos_inc, kTe_hot, kTe_warm, gamma_hot,
                  gamma_warm, r_hot, r_warm, log_rout, fcol, h_max, z)
    
    ragn.set_units('counts')
    ragn.set_flux()
    
    #Ltot = ragn.totSpec_std()
    Ltot = ragn.get_totSED(rel=userel)
    Ld = ragn.get_DiscComponent(rel=userel)
    Lw = ragn.get_WarmComponent(rel=userel)
    Lh = ragn.get_HotComponent(rel=userel)
    nu = ragn.E_obs
    
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    
    ax.loglog(nu, nu**2*Ld, color='red', ls='-.')
    ax.loglog(nu, nu**2*Lw, color='green', ls='-.')
    ax.loglog(nu, nu**2*Lh, color='blue', ls='-.')
    ax.loglog(nu, nu**2*Ltot, color='k')
    
    
    ax.set_ylim(max(nu**2*Ltot)*1e-2, max(nu**2*Ltot)*2)
    plt.show()


if __name__ == '__main__':
    main()



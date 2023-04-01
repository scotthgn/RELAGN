#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 15:55:23 2022

@author: Scott Hagen
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import sys
mdir = os.path.abspath(__file__)
mdir = mdir.replace('/tests/relagn_spec.py', '/src/python_version')
sys.path.append(mdir)

from relagn import relqso


def main():
    #Test pars
    M=1e8
    dist=1
    log_mdot=-1
    a=0.
    cos_inc=0.9
    fcol=1
    z=0
    
    ragn = relqso(M, dist, log_mdot, a, cos_inc, fcol, z)
    print(ragn.r_h, ragn.gamma_h)
    Ltot_n = ragn.totSpec_std()
    Ltot_r = ragn.totSpec_rel()
    nu = ragn.nu_obs
    print(Ltot_r)
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    
    ax.loglog(nu, nu*Ltot_r, color='k')
    ax.loglog(nu, nu*Ltot_n, color='red', ls='dashed')
    
    ax.set_ylim(max(nu*Ltot_n)/1e2, max(nu*Ltot_n)*2)
    plt.show()
    
    
if __name__ == '__main__':
    main()
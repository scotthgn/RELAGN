#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 16:50:31 2022

@author: wljw75
"""
"""
Calculate AGNSED (Kubota & Done 2018) with relativistic corrections at each 
annulus using KYCONV (Dovciak, Karas & Yaqoob 2004)

For the Comptonised region of the disc we use the pyNTHCOMP routine, adapted
from the XSPEC model NTHCOMP (Zdziarski, Johnson & Magdziarz, 1996; Zycki,
Done & Smith, 1999) by Thomas et al. 2016
"""

import numpy as np
from scipy.integrate import quad
import xspec
import warnings
import os
import matplotlib.pyplot as plt

import astropy.constants as const
import astropy.units as u

from pyNTHCOMP import donthcomp


xspec.Xset.chatter = 0 #Stop xspec printing every time it gets called...
#Stop all the run-time warnings (we know why they happen - doesn't affect the output!)
warnings.filterwarnings('ignore') 

"""
Constants
"""
G = (const.G * const.M_sun).value #Grav const in units m^3 s^-1 Msol^-1
c = const.c.value #speed of light, m/s
h = const.h.value #planck constant, Js
sigma_T = const.sigma_T.value #Thompson cross section, m^2
sigma_sb = const.sigma_sb.value #Stefan-Boltzmann const, W m^-2 K^-4
m_p = const.m_p.value #Proton mass, kg
k_B = const.k_B.value #Boltzmann const, J K^-1



"""
Black-Body profile
"""

def do_black_body(T, nu):
    pre_fac = (2 * h * nu**3)/(c**2)
    exp_fac = np.exp((h * nu)/(k_B * T)) - 1
    Bnu = pre_fac / exp_fac

    return np.pi * Bnu


"""
The relagnsed class
"""

class kyagnsed:
    
    Emin = 1e-4
    Emax = 1e4
    numE = 500
    mu = 0.55 #mean particle mass - fixed at solar abundances
    A = 0.3 #Disc albedo = fixed at 0.3 for now
    
    dr_dex = 10 #grid spacing - N points per decade
    
    as_cgs = False #This flags whether to use SI or cgs - False => SI
    as_counts = False #Flags wheter to return as photon counts
    as_flux = False #This flags whether to return luminosity or flux
    
    def __init__(self,
                 M,
                 dist,
                 log_mdot,
                 a,
                 cos_inc,
                 kTe_hot,
                 kTe_warm,
                 gamma_hot,
                 gamma_warm,
                 r_hot,
                 r_warm,
                 log_rout,
                 fcol,
                 h_max,
                 reprocess,
                 z):
        
        """
        Initiates relagnsed object
        
        Parameters
        ----------
        M : float
            Black hole mass - units : Msol
        dist : float
            Co-Moving Distance - units : Mpc
        log_mdot : float
            log mass accretion rate - units : Eddington
        a : float
            Dimensionless Black Hole spin - +ve for prograde rotation,
            -ve for retrograde 
        cos_inc : float
            cos inclination angle
        kTe_hot : float
            Electron temp for hot Compton region - units : keV
        kTe_warm : float
            Electron temp for warm Compton region - units : keV
        gamma_hot : float
            Spectral index for hot Compton region
        gamma_warm : float
            Spectral index for warm Compton region
        r_hot : float
            Outer radius of hot Compton region - units : Rg
        r_warm : float
            Outer radius of warm Compton region - units : Rg
        log_rout : float
            log of outer disc radius - units : Rg
            WARNING! Max value is r_out = 1e3 - set by kyconv
        fcol : float
            Colour temperature correction as described in Done et al. (2012)
            If -ve then follows equation 1 and 2 in Done et al. (2012).
            If +ve then assumes this to be constant correction over entire disc region
        h_max : float
            Scale height of hot Compton region - units : Rg
        reprocess : int
            IF 1 then include re-processing
            IF 0 then dont
        z : float
            Redshift
        """
        
        #Read params
        self.M = M
        self.D, self.d = dist, (dist * u.Mpc).to(u.cm).value
        self.mdot = 10**(log_mdot)
        self.a = np.float64(a)
        self.inc = np.arccos(cos_inc)
        self.kTe_h = kTe_hot
        self.kTe_w = kTe_warm
        self.gamma_h = gamma_hot
        self.gamma_w = gamma_warm
        self.r_h = r_hot
        self.r_w = r_warm
        self.r_out = 10**(log_rout)
        self.fcol = fcol
        self.hmax = h_max 
        self.z = z
        
        self.cosinc = cos_inc
        
        
        #Performing checks
        self._check_spin()
        self._check_inc()
        
        #Calculating disc params 
        self._calc_risco()
        self._calc_r_selfGravity()
        self._calc_Ledd()
        self._calc_efficiency()
        
        
        if log_rout < 0:
            self.r_out = self.r_sg #setting to self gravity if log_rout < 0
        
        if r_warm == -1:
            self.r_w = self.risco
        
        if r_hot == -1:
            self.r_h = self.risco
            reprocess = 0 #No reprocessing if no corona...
        
        self._check_rw()
        self._check_risco()
        self._check_hmax()
        
        #physical conversion factors
        self.Mdot_edd = self.L_edd/(self.eta * c**2)
        self.Rg = (G * self.M)/(c**2)
        self._calc_Dl()
        
        
        #Energy/frequency grid
        self.Egrid = np.geomspace(self.Emin, self.Emax, self.numE)
        self.nu_grid = (self.Egrid * u.keV).to(u.Hz,
                                equivalencies=u.spectral()).value
        self.nu_obs = self.nu_grid/(1 + self.z) #Observers frame
        self.E_obs = self.Egrid/(1 + self.z)
        
        
        #Creating radal grid over disc and warm compton regions
        #using spacing of dr_dex
        self.dlog_r = 1/self.dr_dex
        self.logr_ad_bins = self._make_rbins(np.log10(self.r_w), np.log10(self.r_out))
        self.logr_wc_bins = self._make_rbins(np.log10(self.r_h), np.log10(self.r_w))
        self.logr_hc_bins = self._make_rbins(np.log10(self.risco), np.log10(self.r_h))
        
        #calculating coronal luminosity
        self.reprocess = 0
        self.Lx = self.hotCorona_lumin()
        self.reprocess = reprocess
        
        xspec.AllData.dummyrsp(self.Emin, self.Emax, self.numE)
        
    
    
    
    """
    Performing checks on certain parameters. To ensure that we are within both
    physical limits (i.e -0.998 <= a <= 0.998) AND that we dont wonder off
    acceptable values of kyconv/kerrbb (i.e 3 <= inc <= 85 deg)
    """
    
    def _check_spin(self):
        if self.a >= -0.998 and self.a <= 0.998:
            pass
        else:
            raise ValueError('Spin ' + str(self.a) + ' not physical! \n'
                             'Must be within: -0.998 <= a_star <= 0.998')
    
    
    def _check_inc(self):
        if self.cosinc <= 0.98 and self.cosinc >= 0.09:
            pass
        else:
            raise ValueError('Inclination out of bounds - will not work with kyconv or kerrbb! \n'
                             'Require: 0.09 <= cos(inc) <= 0.98 \n'
                             'Translates to: 11.5 <= inc <= 85 deg')
    
    def _check_rw(self):
        if self.r_w >= self.r_h:
            pass
        else:
            print('WARNING r_warm < r_hot ---- Setting r_warm = r_hot')
            self.r_w = self.r_h
        
        if self.r_w <= self.r_out:
            pass
        else:
            print('WARNING r_warm > r_out ----- Setting r_warm = r_out')
            self.r_w = self.r_out
    
    
    def _check_risco(self):
        if self.r_h >= self.risco:
            pass
        else:
            print('WARNING r_hot < r_isco ----- Setting r_hot = r_isco')
            self.r_h = self.risco
        
        if self.r_w >= self.risco:
            pass
        else:
            print('WARNING! r_warm < r_isco ----- Settin r_warm = r_isco')
            self.r_w = self.risco
    
    
    def _check_hmax(self):
        if self.hmax <= self.r_h:
            pass
        else:
            print('WARNING! hmax > r_h ------- Setting hmax = r_h')
            self.hmax = self.r_h
    
    
    """
    Section for dealing with units. Essentially just methods to change the unit
    flag - and then methods to convert calculated spectrum to desired units.
    Also includes method for changing the energy grid to something other
    than the default
    """
          
    def set_cgs(self):
        """
        Changes output spectra to cgs units
        All outputs will then be in flux of ergs/s/cm^2/Hz
        
        If SI units then outputs in luminosity: W/Hz
        """
        self.as_cgs = True
        self.as_counts = False
    
    def set_counts(self):
        """
        Changes output spectra to photons/s/keV

        """
        self.as_cgs = False
        self.as_counts = True
    
    def set_flux(self):
        """
        Changes output spectra from luminosity to flux

        """
        self.as_flux = True
    
    
    def _to_cgs(self, Lnus):
        """
        Converts input from SI luminosity (W/Hz) to cgs ergs/s

        """
        return Lnus * 1e7
    
    def _to_counts(self, Lnus):
        """
        Converts input from SI luminosity (W/Hz) to photons/s/keV

        """

        flxs = (Lnus * u.W/u.Hz).to(u.keV/u.s/u.keV,
                                                equivalencies=u.spectral()).value
        phs = flxs/self.Egrid
        return phs
        
        
    
    def _to_flux(self, Lnus):
        """
        Converts luminosity to flux. Either per cm^2 or per m^2, depending
        on if SI or cgs

        """
        if self.as_cgs == True or self.as_counts == True:
            dist = self.dl
        else:
            dist = self.dl/100
        
        return Lnus/(4*np.pi*dist**2)
    
    
    def new_ear(self, ear):
        """
        Defines new energy grid if necessary

        Parameters
        ----------
        ear : 1D-array
            New energy grid - units : keV.

        """
        self.Egrid = ear         
        self.nu_grid = (self.Egrid * u.keV).to(u.Hz,
                                equivalencies=u.spectral()).value
        self.nu_obs = self.nu_grid/(1 + self.z) #Observers frame
        self.E_obs = self.Egrid/(1 + self.z)
        
        self.Emin = min(self.Egrid)
        self.Emax = max(self.Egrid)
        self.numE = len(self.Egrid)
        xspec.AllData.dummyrsp(self.Emin, self.Emax, self.numE)
        
        
    
    
    
    """
    Calculating disc properties
    i.e r_isco, L_edd, NT temp, etc
    """
    
    def _calc_Ledd(self):
        """
        Caclulate eddington Luminosity

        """
        #Ledd_c = (4 * np.pi * G * self.M * (2*self.mu) * m_p * c)/sigma_T
        Ledd = 1.39e31 * self.M 
        #print(Ledd, Ledd_c)
        self.L_edd = Ledd
    
    
    def _calc_risco(self):
        """
        Calculating innermost stable circular orbit for a spinning
        black hole. Follows Page and Thorne (1974). Note, can also be reffered
        to as r_ms, for marginally stable orbit
        
        return r_isco as property - so will be called in __init__

        """
        Z1 = 1 + (1 - self.a**2)**(1/3) * (
            (1 + self.a)**(1/3) + (1 - self.a)**(1/3))
        Z2 = np.sqrt(3 * self.a**2 + Z1**2)

        self.risco = 3 + Z2 - np.sign(self.a) * np.sqrt(
            (3 - Z1) * (3 + Z1 + 2*Z2))
    
    
    def _calc_r_selfGravity(self):
        """
        Calcultes the self gravity radius according to Laor & Netzer 1989
        
        NOTE: Assuming that \alpha=0.1 - in future should figure out how to
        constrain this properly!!!

        """
        alpha = 0.1 #assuming turbulence NOT comparable to sound speed
        #See Laor & Netzer 1989 for more details on constraining this parameter
        m9 = self.M/1e9
        self.r_sg = 2150 * m9**(-2/9) * self.mdot**(4/9) * alpha**(2/9)
    
    
    def _calc_efficiency(self):
        """
        Calculates the accretion efficiency eta, s.t L_bol = eta Mdot c^2
        Using the GR case, where eta = 1 - sqrt(1 - 2/(3 r_isco)) 
            Taken from: The Physcis and Evolution of Active Galactic Nuceli,
            H. Netzer, 2013, p.38
        
        Note to self!: When I derive this in Newtonian limit I get
        eta = 1/(2 r_isco). Not entirely sure how to derive the GR version.
        Should ask Chris at next meeting!!!!

        """
        
        self.eta = 1 - np.sqrt(1 - 2/(3*self.risco))
    
    
    def _calc_NTparams(self, r):
        """
        Calculates the Novikov-Thorne relativistic factors.
        see Active Galactic Nuclei, J. H. Krolik, p.151-154
        and Page & Thorne (1974)

        """
        y = np.sqrt(r)
        y_isc = np.sqrt(self.risco)
        y1 = 2 * np.cos((1/3) * np.arccos(self.a) - (np.pi/3))
        y2 = 2 * np.cos((1/3) * np.arccos(self.a) + (np.pi/3))
        y3 = -2 * np.cos((1/3) * np.arccos(self.a))

        
        B = 1 - (3/r) + ((2 * self.a)/(r**(3/2)))
        
        C1 = 1 - (y_isc/y) - ((3 * self.a)/(2 * y)) * np.log(y/y_isc)
        
        C2 = ((3 * (y1 - self.a)**2)/(y*y1 * (y1 - y2) * (y1 - y3))) * np.log(
            (y - y1)/(y_isc - y1))
        C2 += ((3 * (y2 - self.a)**2)/(y*y2 * (y2 - y1) * (y2 - y3))) * np.log(
            (y - y2)/(y_isc - y2))
        C2 += ((3 * (y3 - self.a)**2)/(y*y3 * (y3 - y1) * (y3 - y2))) * np.log(
            (y - y3)/(y_isc - y3))
        
        C = C1 - C2
        
        return C/B
        
        
    
    
    def calc_Tnt(self, r):
        """
        Calculates Novikov-Thorne disc temperature^4 at radius r. 
 
        """
        Rt = self._calc_NTparams(r)
        const_fac = (3 * G * self.M * self.mdot * self.Mdot_edd)/(
            8 * np.pi * sigma_sb * (r * self.Rg)**3)
        
        T4 = const_fac * Rt

        return T4
    
    
    def calc_Trep(self, r):
        """
        Calculates re-processed temperature

        """
        R = r * self.Rg
        H = self.hmax * self.Rg
        
        Frep = (0.5 * self.Lx)/(4*np.pi * (R**2 + H**2))
        Frep *= H/np.sqrt(R**2 + H**2)
        Frep *= (1 - self.A) 
        
        T4rep = Frep/sigma_sb
        return T4rep
    
    
    def calc_Ttot(self, r):
        """
        Caclculates total temperature
        Depends on re-processing flag!
        """
        if self.reprocess == 1:
            T4tot = self.calc_Tnt(r) + self.calc_Trep(r)
        else:
            T4tot = self.calc_Tnt(r)
        
        return T4tot
        
    
    def calc_fcol(self, Tm):
        """
        Calculates colour temperature correction following Eqn. (1) and 
        Eqn. (2) in Done et al. (2012)

        Parameters
        ----------
        Tm : float
            Max temperature at annulus (ie Ttot(r)) - units : K.

        Returns
        -------
        fcol_d : float
            colour temperature correction at T

        """
        if Tm > 1e5:
            #For this region follows Eqn. (1)
            Tm_j = k_B * Tm
            Tm_keV = (Tm_j * u.J).to(u.keV).value #convert to units consitant with equation
            
            fcol_d = (72/Tm_keV)**(1/9)
        
        elif Tm < 1e5 and Tm > 3e4:
            #For this region follows Eqn. (2)
            fcol_d = (Tm/(3e4))**(0.82)
        
        else:
            fcol_d = 1
        
        return fcol_d
    
    
    
    def _calc_Dl(self):
        """
        Calculates luminosity distance to source
        """
        
        self.Dl = self.D * (1+self.z) #In Mpc
        self.dl = self.d * (1+self.z) #In cm
    
    
    
    
    """
    Creating the radial bins to use when calculating each model component
    """
    def _make_rbins(self, logr_in, logr_out):
        """
        Creates an array of radial bin edges, with spacing defined by dr_dex
        Calculates the bin edges from r_out and down to r_in. IF the bin
        between r_in and r_in+dr is less than dlog_r defined by dr_dex, then
        we simply create a slightly wider bin at this point to accomodate
        for the difference

        Parameters
        ----------
        logr_in : float
            Inner radius of model section - units : Rg.
        logr_out : float
            Outer radius of model section - units : Rg.

        Returns
        -------
        logr_bins : 1D-array
            Radial bin edges for section - units : Rg.

        """
        i = logr_out
        logr_bins = np.array([np.float64(logr_out)]) 
        while i > logr_in:
            r_next_edge = i - self.dlog_r
            logr_bins = np.insert(logr_bins, 0, r_next_edge)
            i = r_next_edge

       
        if logr_bins[0] != logr_in:
            if logr_bins[0] < logr_in:
                if len(logr_bins) > 1:
                    logr_bins = np.delete(logr_bins, 0)
                    logr_bins[0] = logr_in
                else:
                    logr_bins[0] = logr_in
            else:
                logr_bins[0] = logr_in
        
        return logr_bins
        
    
    
        
    
    
    
    
    
    """
    Section to calculate annulus spectra
    """
    
    def disc_annuli(self, r, dr):
        """
        Calculates disc spectrum for annulus at position r with width dr. Note
        that r is taken to be the center of the bin!

        Parameters
        ----------
        r : float
            Inner radius of annulus - units : Rg.
        dr : float
            Width of annulus - units : Rg.
        
        Returns
        -------
        Lnu_ann : 1D-array
            Disc black-body at annulus - units : W/Hz

        """
        T4_ann = self.calc_Ttot(r)
        Tann = T4_ann**(1/4)
        if self.fcol < 0:
            fcol_r = self.calc_fcol(Tann)
        else:
            fcol_r = self.fcol
        
        Tann *= fcol_r

        #Extend nu/energy grid to avoid inaccuracies at large radii
        #Only do this if internal Emin > 1e-4 keV
        if self.Emin < 1e-4:
            nu_ext = np.insert(self.nu_grid, 0, np.geomspace(1e12, min(self.nu_grid)-100, 50))
            nu_use = nu_ext
        
        else:
            nu_use = self.nu_grid
        
        B = do_black_body(Tann, nu_use)
        
        #Annulus luminosity - normalised as black body
        #NOTE! kyconv integrates over the annulus within the code, 
        #hence why we are NOT multiplying with 4pi r dr - as this would 
        #effectively overestimate the actual normalisation post convolution.
        #The final result is still the same as this - it's just done somewhere else...
        #kyconv also deal with the inclination - so no cos(inc) correction factor needs to be applied
        norm = sigma_sb * (Tann/fcol_r)**4 * self.Rg**2

        radiance = np.trapz(B, nu_use)
        if radiance == 0:
            Lnu_ann = np.zeros(len(nu_use))
        else:
            Lnu_ann = norm * (B/radiance)
        
        #Re-casting onto original grid in order to be consistent with kyconv
        #i.e, just slicing away entries below Emin/nu_grid_min
        if self.Emin < 1e-4:
            Lnu_ann = Lnu_ann[50:]
            
        return Lnu_ann
        
    
    
    def warmComp_annuli(self, r, dr):
        """
        Calculates comptonised spectrum for annulus at r with width dr. Note,
        r taken to be in center of bin!
        Uses pyNTHCOMP

        Parameters
        ----------
        r : float
            Inner radius of annulus - units : Rg.
        dr : float
            Width of annulus - units : Rg.
        
        Returns
        -------
        Lnu_ann : 1D-array
            Warm Comptonised spectrum at annulus - units : W/Hz
        """
        
        T4_ann = self.calc_Ttot(r)
        Tann = T4_ann**(1/4)
        
        kTann = k_B * Tann
        kTann = (kTann * u.J).to(u.keV).value #converting T to keV for nthcomp
        
        ph_nth = donthcomp(self.Egrid, [self.gamma_w, self.kTe_w,
                                        kTann, 1, 0])
        ph_nth = (ph_nth * u.W/u.keV).to(u.W/u.Hz, 
                                            equivalencies=u.spectral()).value
        
        norm = sigma_sb * (Tann**4) * self.Rg**2
        radiance = np.trapz(ph_nth, self.nu_grid)
        if radiance == 0:
            Lnu_ann = np.zeros(len(self.nu_grid))
        else:
            Lnu_ann = norm * (ph_nth/radiance)# * self.cosinc/0.5
        
        return Lnu_ann
    
    
    
    
    
    
    
    """
    Section for calculating total spectrum for each disc component
    So spectrum for disc component and warm compton component
    Including methods for both with and without kyconv
    """
    
    
    def do_relDiscSpec(self):
        """
        Calculates contribution from entire disc section - for relativistic 
        case

        """    
        for i in range(len(self.logr_ad_bins) - 1):
            dr_bin = 10**self.logr_ad_bins[i+1] - 10**self.logr_ad_bins[i] #width in lin space
            rmid = 10**(self.logr_ad_bins[i] + self.dlog_r/2) #geometric center of bin
            Lnu_ann = self.disc_annuli(rmid, dr_bin)

                
                #Defining XSPEC model in order to convolve with kyconv
            def disc_ann(es, params, flx):
                Els = np.array(es[:-1])
                Ers = np.array(es[1:])
                
                dEs = Ers - Els
                Emids = Els + dEs/2
                
                #Multiply by 4 because kyconv only does pi r dr
                #We want 4 pi r dr
                fluxs = (Lnu_ann*4 * u.W/u.Hz).to(u.keV/u.s/u.keV,
                                                equivalencies=u.spectral()).value
                
            
                fluxs = fluxs/(4*np.pi * self.dl**2)
    
                phs = fluxs/Emids
                for j in range(len(es) - 1):
                    flx[j] = dEs[j]*phs[j]
                        
                        
            parinfo_ad = ('pn "" 1 0.1 0 0 2 2',) #just a dummy parameter to make xspec happy
            xspec.AllModels.addPyMod(disc_ann, parinfo_ad, 'add')
            
            #Convolving annulus with kyconv
            #kyconv params are:
                #a, spin
                #inc - deg
                #r_in - Rg
                #ms - set to 0 for integration to start at r_in
                #r_out - Rg
                #alph - Emissivity power law - set to 0 as this already hardwired in own code
                #beta - Same as alpha
                #rb - break radius for emissivity - dummy parameter in this scenario
                #z - redshift, assume 0 for now
                #ntable - set to 0 for isotropic emission
                #nE - Energy resolution - set to same as rest of model
                #norm - set to -1 for NO renormalisation - maintain physical units
            if self.logr_ad_bins[i] <= 3 and self.logr_ad_bins[i+1] <= 3:
                r_in = 10**self.logr_ad_bins[i]
                r_out = 10**self.logr_ad_bins[i+1]
                r_br = rmid
                
                kyparams = (self.a, np.rad2deg(self.inc), r_in, 0, r_out, 0, 0,
                         r_br, 0, 0, self.numE, -1)
                
                xspec.Model('kyconv*disc_ann', setPars=kyparams)
                #xspec.Model('disc_ann')
                xspec.Plot.device = '/null'
                
                xspec.Plot('model')
                Es = np.array(xspec.Plot.x())
                phs_r = np.array(xspec.Plot.model())
                
                L_kev = phs_r * Es * 4 * np.pi * self.dl**2
                Lnu_r = (L_kev * u.keV/u.s/u.keV).to(u.W/u.Hz, equivalencies=u.spectral()).value
                
                xspec.AllModels.clear()
                
            else:
                #If out of bounds for tables do non-relativistic
                #This is fine - as beyon 1000Rg Gr effects tiny!
                #Since no kyconv - we now apply the normalisation here instead
                Lnu_r = Lnu_ann * 2*np.pi * 2 * rmid * dr_bin   * self.cosinc/0.5
            

            if i == 0:
                Lnu_all = Lnu_r
            else:
                Lnu_all = np.column_stack((Lnu_all, Lnu_r))
                
        if np.shape(Lnu_all) != np.shape(self.Egrid):
            Lnu_tot = np.sum(Lnu_all, axis=-1)
        else:
            Lnu_tot = Lnu_all
        
        self.Lnu_disc_rel = Lnu_tot  #* self.cosinc/0.5
        return Lnu_tot
    
    
    def do_stdDiscSpec(self):
        """
        Calculates contribution from entire disc section - for non-relativisitc
        case. Usefull for comparison...

        """
        for i in range(len(self.logr_ad_bins) - 1):
            dr_bin = 10**self.logr_ad_bins[i+1] - 10**self.logr_ad_bins[i]
            rmid = 10**(self.logr_ad_bins[i] + self.dlog_r/2)
            
            Lnu_r = self.disc_annuli(rmid, dr_bin) * 2*np.pi * 2 * rmid * dr_bin
            
            if i == 0:
                Lnu_all = Lnu_r
            else:
                Lnu_all = np.column_stack((Lnu_all, Lnu_r))
                
        if np.shape(Lnu_all) != np.shape(self.Egrid):
            Lnu_tot = np.sum(Lnu_all, axis=-1)
        else:
            Lnu_tot = Lnu_all
        
        self.Lnu_disc_norel = Lnu_tot  * self.cosinc/0.5
        return self.Lnu_disc_norel
            
    
    
    def do_relWarmCompSpec(self):
        """
        Calculates contribution from entire warm Compton region - for 
        relativistic case

        """
        for i in range(len(self.logr_wc_bins) - 1):
            dr_bin = 10**self.logr_wc_bins[i+1] - 10**self.logr_wc_bins[i]
            rmid = 10**(self.logr_wc_bins[i] + self.dlog_r/2)
            
            Lnu_ann = self.warmComp_annuli(rmid, dr_bin)
            
            #creating pyxspec model so can do convolution with kyconv
            def warm_ann(es, params, flx):
                Els = np.array(es[:-1])
                Ers = np.array(es[1:])
                
                dEs = Ers - Els
                Emids = Els + dEs/2
                
                fluxs = (Lnu_ann*4 * u.W/u.Hz).to(u.keV/u.s/u.keV,
                                                equivalencies=u.spectral()).value
                
            
                fluxs = fluxs/(4*np.pi * self.dl**2)
                
                phs = fluxs/Emids
                for j in range(len(es) - 1):
                    flx[j] = dEs[j]*phs[j]
            
            parinfo_ad = ('pn "" 1 0.1 0 0 2 2',) #just a dummy parameter to make xspec happy
            xspec.AllModels.addPyMod(warm_ann, parinfo_ad, 'add')
            
            #Convolving annulus with kyconv
            #See do_relDiscSpec for explanation of kyparams
            if self.logr_wc_bins[i] <= 3 and self.logr_wc_bins[i+1] <= 3:
                r_in = 10**self.logr_wc_bins[i]
                r_out = 10**self.logr_wc_bins[i+1]
                r_br = rmid
                kyparams = (self.a, np.rad2deg(self.inc), r_in, 0, r_out, 0, 0,
                            r_br, 0, 0, self.numE, -1)

                xspec.Model('kyconv*warm_ann', setPars=kyparams)
                xspec.Plot.device = '/null'
                
                xspec.Plot('model')
                Es = np.array(xspec.Plot.x())
                phs_r = np.array(xspec.Plot.model())
                
                
                L_kev = phs_r * Es * 4 * np.pi * self.dl**2
                Lnu_r = (L_kev * u.keV/u.s/u.keV).to(u.W/u.Hz, equivalencies=u.spectral()).value
                
                xspec.AllModels.clear()
            
            else:
                #If out of bounds for tables do non-relativistic
                #This is fine - as beyon 1000Rg Gr effects tiny!
                #Since no kyconv - we now apply the normalisation here instead
                Lnu_r = Lnu_ann * 2*np.pi * 2 * rmid * dr_bin   * self.cosinc/0.5
                

            if i == 0:
                Lnu_all = Lnu_r
            else:
                Lnu_all = np.column_stack((Lnu_all, Lnu_r))
        
        if np.shape(Lnu_all) != np.shape(self.Egrid):
            Lnu_tot = np.sum(Lnu_all, axis=-1)
        else:
            Lnu_tot = Lnu_all
            
        self.Lnu_warm_rel = Lnu_tot
        return Lnu_tot
    
    
    def do_stdWarmCompSpec(self):
        """
        Calculates contribution from entire warm Compton region - for 
        non-relativistic case
        
        """
        for i in range(len(self.logr_wc_bins) - 1):
            dr_bin = 10**self.logr_wc_bins[i+1] - 10**self.logr_wc_bins[i]
            rmid = 10**(self.logr_wc_bins[i] + self.dlog_r/2)
            
            Lnu_r = self.warmComp_annuli(rmid, dr_bin) * 4*np.pi*rmid*dr_bin
            
            if i == 0:
                Lnu_all = Lnu_r
            else:
                Lnu_all = np.column_stack((Lnu_all, Lnu_r))
        
        if np.shape(Lnu_all) != np.shape(self.Egrid):
            Lnu_tot = np.sum(Lnu_all, axis=-1)
        else:
            Lnu_tot = Lnu_all
            
        self.Lnu_warm_norel = Lnu_tot * self.cosinc/0.5
        return self.Lnu_warm_norel
        
        
    
    """
    Section for calculating the hot coronal emission
    """
    def seed_tempHot(self):
        """
        Calculated seed photon temperature for the hot compton region.
        Follows xspec model agnsed, from Kubota & Done (2018)
        
        Returns
        -------
        kT_seed : float
            Seed photon temperature for hot compton - units : keV

        """
        T4_edge = self.calc_Ttot(self.r_h) #inner disc T in K
        Tedge = T4_edge**(1/4)
        
        kT_edge = k_B * Tedge #units J
        kT_edge = (kT_edge * u.J).to(u.keV).value
        if self.r_w != self.r_h:
            #If there is a warm compton region then seed photons mostly from
            #here. Will then need to include Compton y-param
            ysb = (self.gamma_w * (4/9))**(-4.5)
            kT_seed = np.exp(ysb) * kT_edge
            #kT_seed = kT_edge
            
        else:
            #If only disc down to r_hot seed temp will be same as inner disc temp
            kT_seed = kT_edge
            if self.fcol < 0:
                fcol_r = self.calc_fcol(Tedge)
            else:
                fcol_r = self.fcol
            
            kT_seed *= fcol_r
        
        return kT_seed
    
    
    def Lseed_hotCorona(self):
        """
        Calculated luminsoty of seed photons emitted at radius r, intercepted
        by corona

        Returns
        -------
        Lseed_tot : float
            Total seed photon luminosity seen by corona - units : W

        """
        logr_all_grid = self._make_rbins(np.log10(self.r_h), np.log10(self.r_out))
        Lseed_tot = 0
        for i in range(len(logr_all_grid) - 1):
            dr = 10**logr_all_grid[i+1] - 10**logr_all_grid[i]
            rmid = 10**(logr_all_grid[i] + 0.5 * self.dlog_r)
            
            if self.hmax <= rmid:
                theta_0 = np.arcsin(self.hmax/rmid)
                cov_frac = theta_0 - 0.5 * np.sin(2*theta_0)
            else:
                cov_frac = 1
            
            T4_ann = self.calc_Ttot(rmid)
            
            Fr = sigma_sb * T4_ann
            Lr = 2 * 2*np.pi*rmid*dr * Fr * cov_frac/np.pi * self.Rg**2


            Lseed_tot += Lr
        
        return Lseed_tot
        
    
    
    def hotCorona_lumin(self):
        """
        Calculates the coronal luminosity - used as normalisaton for the 
        hot compton spectral component.
        
        Calculated as Lhot = Ldiss + Lseed
        
        where Ldiss is the energy dissipated from the accretion flow, and
        Lseed is the seed photon luminosity intercpted by the corona

        """
        
        Ldiss, err = quad(lambda rc: 2*sigma_sb*self.calc_Tnt(rc) * 2*np.pi*rc * self.Rg**2,
                     self.risco, self.r_h)

        Lseed = self.Lseed_hotCorona()
        
        Lhot = Ldiss + Lseed
        return Lhot
    
    
    def hotComp_annuli(self, r, dr):
        """
        Calculates spectrum from radial slice of hot comp region
        Neccessary in order to apply kyconv correctly!

        """
        kT_ann = self.seed_tempHot()
        T4ann = self.calc_Tnt(r)
        Ldiss_ann = sigma_sb * T4ann * self.Rg**2
        #For Lseed need to /4pi r dr so that normalise correctly in kyconv
        #Also need to /N_bins - assuming the corona has no radial dependence in its emissivity
        Lseed_ann = self.Lseed_hotCorona()/(4*np.pi*r*dr * (len(self.logr_hc_bins) - 1))
        Ltot_ann = Ldiss_ann + Lseed_ann
        
        if kT_ann < self.kTe_h:
            ph_ann = donthcomp(self.Egrid, [self.gamma_h, self.kTe_h, kT_ann, 1, 0])
            ph_ann = (ph_ann * u.W/u.keV).to(u.W/u.Hz, 
                                            equivalencies=u.spectral()).value
            Lnu_ann = Ltot_ann * (ph_ann/np.trapz(ph_ann, self.nu_grid))
            
        else:
            Lnu_ann = np.zeros(len(self.Egrid))
            
        return Lnu_ann
        
    
    
    def do_relHotCompSpec(self):
        """
        Calculates spectrum of hot compton region - with relativity!

        """
        for i in range(len(self.logr_hc_bins) - 1):
            dr_bin = 10**self.logr_hc_bins[i+1] - 10**self.logr_hc_bins[i]
            rmid = 10**(self.logr_hc_bins[i] + self.dlog_r/2)
            
            Lnu_ann = self.hotComp_annuli(rmid, dr_bin)
            
            #creating pyxspec model so can do convolution with kyconv
            def hot_ann(es, params, flx):
                Els = np.array(es[:-1])
                Ers = np.array(es[1:])
                
                dEs = Ers - Els
                Emids = Els + dEs/2
                
                fluxs = (Lnu_ann*4 * u.W/u.Hz).to(u.keV/u.s/u.keV,
                                                equivalencies=u.spectral()).value
                
            
                fluxs = fluxs/(4*np.pi * self.dl**2)
                
                phs = fluxs/Emids
                
                for j in range(len(es) - 1):
                    flx[j] = dEs[j]*phs[j]
            
            parinfo_ad = ('pn "" 1 0.1 0 0 2 2',) #just a dummy parameter to make xspec happy
            xspec.AllModels.addPyMod(hot_ann, parinfo_ad, 'add')
            
            if self.logr_hc_bins[i] <= 3 and self.logr_hc_bins[i+1] <= 3:
                r_in = 10**self.logr_hc_bins[i]
                r_out = 10**self.logr_hc_bins[i+1]
                r_br = rmid
                
                #Convolving annulus with kyconv
                kyparams = (self.a, np.rad2deg(self.inc), r_in, 0, r_out, 0, 0,
                            r_br, 0, 0, self.numE, -1)
                
                xspec.Model('kyconv*hot_ann', setPars=kyparams)
                xspec.Plot.device = '/null'
                
                xspec.Plot('model')
                Es = np.array(xspec.Plot.x())
                phs_r = np.array(xspec.Plot.model())
                
                
                L_kev = phs_r * Es * 4 * np.pi * self.dl**2
                Lnu_r = (L_kev * u.keV/u.s/u.keV).to(u.W/u.Hz, equivalencies=u.spectral()).value
                
                xspec.AllModels.clear()
            
            else:
                #If out of bounds for tables do non-relativistic
                #This is fine - as beyon 1000Rg Gr effects tiny!
                #Since no kyconv - we now apply the normalisation here instead
                Lnu_r = Lnu_ann * 2*np.pi * 2 * rmid * dr_bin   * self.cosinc/0.5

            if i == 0:
                Lnu_all = Lnu_r
            else:
                Lnu_all = np.column_stack((Lnu_all, Lnu_r))
        
        
        if np.shape(Lnu_all) != np.shape(self.Egrid):
            Lnu_tot = np.sum(Lnu_all, axis=-1)
        else:
            Lnu_tot = Lnu_all
            
        self.Lnu_hot_rel = Lnu_tot
        return Lnu_tot
        
    
    def do_stdHotCompSpec(self):
        """
        Calculates spectrum of hot comptonised region - no relativity

        """
        kTseed = self.seed_tempHot()
        Lum = self.hotCorona_lumin()
        
        if kTseed < self.kTe_h:
            ph_hot = donthcomp(self.Egrid, [self.gamma_h, self.kTe_h, kTseed, 1, 0])
            ph_hot = (ph_hot * u.W/u.keV).to(u.W/u.Hz, 
                                            equivalencies=u.spectral()).value
        
            self.Lnu_hot_norel = Lum * (ph_hot/np.trapz(ph_hot, self.nu_grid))
        
        else:
            self.Lnu_hot_norel = np.zeros(len(self.Egrid))
            
        return self.Lnu_hot_norel 
        
            
    
    
    """
    Return the total sepctrum
    Sum of all components
    
    Unless de-buggin these functions should be used for all things spectral!
    Since they also deal with awkwardness when you set strange limits
    """
    def totSpec_rel(self):
        """
        Total spectrum from the relativistic components
        
        """
        
        #Checking if components already exist in order to avoid multiple calculations
        if hasattr(self, 'Lnu_disc_rel'):
            Ld = self.Lnu_disc_rel
        elif self.r_w != self.r_out and len(self.logr_ad_bins) != 1:
            Ld = self.do_relDiscSpec()
        else:
            Ld = np.zeros(len(self.nu_grid))
        
        if hasattr(self, 'Lnu_warm_rel'):
            Lw = self.Lnu_warm_rel
        elif self.r_w != self.r_h and len(self.logr_wc_bins) != 1:
            Lw = self.do_relWarmCompSpec()
        else:
            Lw = np.zeros(len(self.nu_grid))
        
        if hasattr(self, 'Lnu_hot'):
            Lh = self.Lnu_hot
        elif self.r_h != self.risco and len(self.logr_hc_bins) != 1:
            Lh = self.do_relHotCompSpec()
        else:
            Lh = np.zeros(len(self.nu_grid))
        
        Lnu_tot_rel = Ld + Lw + Lh
 
        
        #Handling unit conversions if necessary
        if self.as_cgs == True:
            Lnu_tot_rel = self._to_cgs(Lnu_tot_rel)
            Ld = self._to_cgs(Ld)
            Lw = self._to_cgs(Lw)
            Lh = self._to_cgs(Lh)
        
        elif self.as_counts == True:
            Lnu_tot_rel = self._to_counts(Lnu_tot_rel)
            Ld = self._to_counts(Ld)
            Lw = self._to_counts(Lw)
            Lh = self._to_counts(Lh)
        
        
        if self.as_flux == True:
            Lnu_tot_rel = self._to_flux(Lnu_tot_rel)
            Ld = self._to_flux(Ld)
            Lw = self._to_flux(Lw)
            Lh = self._to_flux(Lh)
            
        
        self.Lnu_tot_rel = Lnu_tot_rel
        self.Lnu_disc_rel = Ld
        self.Lnu_warm_rel = Lw
        self.Lnu_hot_rel = Lh
        return self.Lnu_tot_rel
    
    
    def totSpec_std(self):
        """
        Total spectrum from the non-relativistic components

        """
        #Checking if components already exist in order to avoid multiple calculations
        if hasattr(self, 'Lnu_disc_norel'):
            Ld = self.Lnu_disc_norel
        elif self.r_w != self.r_out and len(self.logr_ad_bins) != 1:
            Ld = self.do_stdDiscSpec()
        else:
            Ld = np.zeros(len(self.nu_grid))
        
        if hasattr(self, 'Lnu_warm_norel'):
            Lw = self.Lnu_warm_norel
        elif self.r_w != self.r_h and len(self.logr_wc_bins) != 1:
            Lw = self.do_stdWarmCompSpec()
        else:
            Lw = np.zeros(len(self.nu_grid))
        
        if hasattr(self, 'Lnu_hot_norel'):
            Lh = self.Lnu_hot_norel
        elif self.r_h != self.risco and len(self.logr_hc_bins) != 1:
            Lh = self.do_stdHotCompSpec()
        else:
            Lh = np.zeros(len(self.nu_grid))
        
        Lnu_tot_norel = Ld + Lw + Lh
        
        #Handling unit conversions if necessary
        if self.as_cgs == True:
            Lnu_tot_norel = self._to_cgs(Lnu_tot_norel)
            Ld = self._to_cgs(Ld)
            Lw = self._to_cgs(Lw)
            Lh = self._to_cgs(Lh)
        
        elif self.as_counts == True:
            Lnu_tot_norel = self._to_counts(Lnu_tot_norel)
            Ld = self._to_counts(Ld)
            Lw = self._to_counts(Lw)
            Lh = self._to_counts(Lh)
            
        
        if self.as_flux == True:
            Lnu_tot_norel = self._to_flux(Lnu_tot_norel)
            Ld = self._to_flux(Ld)
            Lw = self._to_flux(Lw)
            Lh = self._to_flux(Lh)
            
        
        self.Lnu_tot_norel = Lnu_tot_norel
        self.Lnu_disc_norel = Ld
        self.Lnu_warm_norel = Lw
        self.Lnu_hot_norel = Lh
        return self.Lnu_tot_norel
    



    
    
    
    
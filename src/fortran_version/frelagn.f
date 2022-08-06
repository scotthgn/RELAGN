c      Fortran version of RELAGN
c      Ref. Hagen and Done (in prep.)
c
c      Calculates AGNSED (Kubota & Done 2018) with relativistic
c      corrections applied at each annulus using KYCONV
c      (Dovciak, Karas & Yaqoob 2004)
c
c      The model consists of three regions: a standard outer disc, a warm
c      Comptonization region where the disc has not thermalised properly,
c      and a spherical inner corona modelled as a hot Comptonization
c      region. All the energetics are calculated in the rest frame of the
c      black hole, before the relativistic transfer functions are applied
c      The Comptonization is modelled using NTHCOMP (Zdziarski, Johnson &
c      Magdziarz 1996; Zycki, Done & Smith 1999)
c     
c      Note to self: start all fortran lines in column 7!!!!!!
c      Goodnes knows why!? I miss python :(
      
       subroutine relagn(ear,ne,param,ifl,photar,photer)
       implicit none

       integer npars
       parameter(npars=15)
      
       integer ne, ifl
       real ear(0:ne), photar(ne), param(npars), photer(ne)
       real oldpar(npars)

       integer Nnew !Defining own energy grid for model calc.
       parameter(Nnew=2000) !Number of bin edges
       real enew(0:Nnew), ph(Nnew) !ph is photon array on new grid
       real dloge, newemin, newemax

       real fstart(ne), fend(ne)
       integer istart(ne), iend(ne)

       logical parchange, echange
       save oldpar

       integer i, n !Iteration indeces

c      ear: Energy array
c      ne: Size of energy array
c      param: Parameter value array
c      ifl: Spectrum number
c      photar: flux array (size ne-1)
c
c      param(1):  Mass, Msol
c      param(2):  Distance, Mpc
c      param(3):  log mdot, Mass accretion rate,  L/Ledd
c      param(4):  astar, BH spin
c      param(5):  cos inc, Inclination
c      param(6):  kTe_hot, hot corona temperature, keV
c      param(7):  kTe_warm, warm compton temperature, keV
c      param(8):  gamma_hot, hot corona photon index
c      param(9):  gamma_warm, warm compton photon index
c      param(10): r_hot, hot corona radius, Rg
c      param(11): r_warm, warm compton outer radius, Rg
c      param(12): log r_out, outer disc radius, Rg
c      param(13): fcol, colour temperature correction - Disc region only!
c      param(14): hmax, scale height of corona, Rg
c      param(15): redshift

c      checking if parameters have changed
       parchange=.false.
       do i=1,npars,1
          if (param(i).ne.oldpar(i)) parchange=.true.
       end do
       
c      Checking if we need to change energy grid
c      The default extends from 1e-4 to 1e3 keV
       if (ear(0).eq.0.0) then
          newemin=ear(1) - (ear(2)-ear(1))/10.0
       else
          newemin=min(1.0e-4, ear(0))
       end if
       newemax=max(1.0e3, ear(ne))

       if ((enew(0).ne.newemin).or.(enew(Nnew).ne.newemax)) then
          echange=.true.
       end if
       
       
c      Calculating new energy grid if necessary
       if (echange) then
          dloge=log10(newemax/newemin)/float(Nnew)
          enew(0) = newemin
          do n=1, Nnew, 1
             enew(n)=10**(log10(enew(0))+dloge*float(n))
          end do
       end if

       call calc_spec(enew, Nnew, param, ifl, ph)

c      Re-bin onto original xspec grid
       call inibin(Nnew, enew, ne, ear, istart, iend, fstart, fend, 0)
       call erebin(Nnew, ph, ne, istart, iend, fstart, fend, photar)

       return
       end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Main subroutine for calculating disc spec, warm Comp spec, and
c      hot Comp spec
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


       subroutine calc_spec(es, nn, param, ifl, ph)

c      Function to calculate the SED. Calculates the spectral
c      contribution from each region (disc, warm comp, hot comp)
c      Temperature found through the Novikov-Thorne equiations
c       
c      The relativistic corrections are applied at each annulus through
c      kyconv

       implicit none
       integer nn, ifl
       real es(0:nn), ph(nn), param(*)
       real ph_ann(nn)
       double precision M, mdot
       double precision rout, rw, rh, rsg, risco
       double precision astar, hmax
       double precision fcol, kTh, kTw
       double precision gammah, gammaw
       double precision cos_inc, dist

       double precision pi, G, c, h, sigma_sb,  mp, kB
       double precision eta, Ledd, Mdot_edd, Rg
       double precision dr_dex, nrd, nrw, nrh
       double precision dlog_rd, dlog_rw, dlog_rh
       double precision rmid, dr
       double precision kkev, kevhz
       double precision dflux, en, dnphot, tr

       double precision calc_risco, calc_rsg
       double precision efficiency, nttemp

       real kpar(12), ntwpar(5), nthpar(5)
       real kpherr(nn), ntwpherr(nn), nthpherr(nn)

       logical return_disc, return_warm, return_hot
       integer i, n
       
c      Constants
       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33 !Grav const. cm^-3 s^-1 Msol^-1
       c = 3.0d10            !Speed of light. cm/s
       h = 6.62d-27          !Plank const. erg/s
       sigma_sb = 5.67d-5    !Stefan-Boltzmann const. erg s^-1 cm^-2 K^-4
       mp = 1.67d-24         !Proton mass. g
       kB = 1.38d-16         !Boltzmann const. erg/K

c      Unit conversion constants
       kkev = 1.6048d7          !K/keV
       kevhz = 2.417965d17      !Hz/keV

       dr_dex = 30              !Radial resolution, Nbins per decade (ish)
      

c      Read parameters
       M = dble(param(1))
       dist = dble(param(2)) !Mpc
       mdot = dble(10**(param(3)))
       astar = dble(param(4))
       cos_inc = dble(param(5))
       kTh = dble(abs(param(6)))
       kTw = dble(abs(param(7)))
       gammah = dble(param(8))
       gammaw = dble(abs(param(9)))
       rw = dble(param(11))
       fcol = dble(param(13))
       hmax = dble(param(14))
       
       
c      Getting accretion attributes
       risco = calc_risco(astar)
       eta = efficiency(risco)
       rsg = calc_rsg(M, mdot)

       Ledd = 1.39d38 * M       !erg/s
       Mdot_edd = Ledd/(eta * c**2) !g/s
       Rg = (G*M)/c**2 !cm
       write(*,*) Mdot_edd * mdot
       dist = dist*1.0d6        !pc
       dist = dist*3.086d18    !cm
       
c      Filling parameter arrays (for those that dont change!!!)
c      Starting with kyconv params
       kpar(1) = astar
       kpar(2) = acos(cos_inc) * (180.0/pi)
       kpar(4) = 0 !set to integrate from r_in (inner edge of annulus)
       kpar(6) = 0 !emissivity power law - assume constant over annulus
       kpar(7) = 0 !same as parameter 6
       kpar(9) = 0 !redshift - dealt with within this script
       kpar(10) = 0 !Limb-darkening. 0=>isotropic
       kpar(11) = nn !Set energy resolution to same as grid
       kpar(12) = -1 !No re-normalisation => want physical!!
c      Note kpar(3)=r_in, kpar(5)=r_outm and kpar(8)=r_br
c      These are set seperately for each annulus!!!!!
       
c      Now warm nthcomp params
       ntwpar(1) = gammaw
       ntwpar(2) = kTw
       ntwpar(4) = 0 !BB assumed from each annulus
       ntwpar(5) = 0 !Redshift - dealt with within this script

c      And finally hot nthomp params
       nthpar(1) = gammah
       nthpar(2) = kTh
       nthpar(4) = 0
       nthpar(5) = 0

       

c      Checking switching parameters
       return_disc=.true.
       return_warm=.true.
       return_hot=.true.
       
       if (param(12).lt.0.0) then !checkig outer radius
          rout = rsg
       else
          rout = dble(10**param(12))
       end if

       if ((param(10).lt.0.0).or.(param(10).lt.risco)) then !Checking inner disc radius
          rh = risco
          return_hot=.false.
       else
          rh = dble(param(10))
       end if

       if ((param(11).lt.0.0).or.(param(11).lt.risco)) then !Checking warm region
          rw = risco
          return_warm=.false.
          return_hot=.false.
       else
          rw = dble(param(11))
       end if

       
       if (param(9).lt.0.0) then
          return_warm=.false.
          return_hot=.false.
       end if

       if (param(7).lt.0.0) then
          return_disc=.false.
          return_hot=.false.
       end if

       if (param(6).lt.0.0) then
          return_disc=.false.
          return_warm=.false.
       end if


       if (rw.lt.rh) then
          call xwrite('r_warm < r_hot!!! Re-setting r_w = r_h', 10)
          rw = rh
          return_warm=.false.
       end if

       if (rw.gt.rout) then
          call xwrite('r_warm > r_out!!! Re-setting r_w = r_out', 10)
          rw = rout
          return_disc=.false.
       end if


c-----------------------------------------------------------------------
c      Section for calculating the disc region
c-----------------------------------------------------------------------
       if (return_disc) then
          nrd = (log10(rout) - log10(rw)) * dr_dex !nr of bins in region
          nrd = ceiling(nrd) !round up to nearest integer value
          dlog_rd = (log10(rout) - log10(rw))/nrd
          do i=1, int(nrd), 1
             rmid = 10**(log10(rw)+float(i-1)*dlog_rd+dlog_rd/2.0)
             dr = 10**(log10(rmid)+dlog_rd/2.0)
             dr = dr - 10**(log10(rmid)-dlog_rd/2.0)

             if ((rmid+dr/2.0 + dr).gt.rout) then !Makes big bin at edge
                dr = rout - 10**(log10(rmid)-dlog_rd/2.0)
             end if

             tr = nttemp(rmid, M, mdot, Mdot_edd, astar, risco) !Temp. K
c            Calculating BB emission over each energy
             do n=1, nn, 1
                en = dble(log10(es(n-1)) + log10(es(n)))
                en = en/2.0 !Evaluate in mid bin
                en = 10**en
                
                if (en.lt.30.0*tr) then
                   dflux = (pi*2.0 * h * (en*kevhz)**3)/(c**2)
                   dflux = dflux * 1/(exp((h*en*kevhz)/(kB*tr)) - 1)
                   dflux = dflux * 4.0 *  Rg**2 !Lum in erg/s/Hz
                   !dflux = dflux * pi*rmid*dr
c                  Note, not applyin pi*r*dr as kyconv deals with this
c                  kyconv also does cos(inc)/0.5

                   dflux = dflux/(4.0*pi*dist**2) !erg/s/cm^2/Hz
                   
                   ph_ann(n) = dflux/(h*kevhz*en) !photons/s/cm^2/Hz
                   ph_ann(n) = ph_ann(n)*kevhz*(es(n)-es(n-1)) !photons/s/cm^2                
                else
                   ph_ann(n) = 0.0d0
                end if         
             end do

c            Now applying kyconv (IF r within the data tables)
             if ((rmid+dr/2.0).lt.1.0d3) then
                kpar(3) = rmid - dr/2.0
                kpar(5) = rmid + dr/2.0
                kpar(8) = rmid
                call kyconv(es, nn, kpar, ifl, ph_ann, kpherr)
             else
                ph_ann = ph_ann * pi*rmid*dr*(cos_inc/0.5)
             end if

c            Adding into main photon array
             do n=1, nn, 1
                if (i.eq.1) then
                   ph(n) = ph_ann(n) 
                else
                   ph(n)=ph(n)+ph_ann(n)
                end if
             end do 
          end do
       end if

          

       end
             






       
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Functions for calculating disc properties
c      e.g risco, eta (efficiency), rsg, disc T, etc.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

       function calc_risco(astar)
c      Calculates the innermost stable circular orbit/marginally stable
c      orbit
       implicit none

       double precision astar, calc_risco
       double precision Z1_1, Z1_2, Z1, Z2
       
       Z1_1 = (1.0 - astar**2.0)**(1.0/3.0)
       Z1_2 = (1.0 + astar)**(1.0/3.0) + (1.0 - astar)**(1.0/3.0)
       Z1 = 1.0 + Z1_1 * Z1_2

       Z2 = sqrt(3.0 * astar**2.0 + Z1**2.0)
    
       if(astar.ge.0.0) then
          calc_risco = 3.0 + Z2 - sqrt((3.0-Z1) * (3.0+Z1+2.0*Z2))
       else
          calc_risco = 3.0 + Z2 + sqrt((3.0-Z1) * (3.0+Z1+2.0*Z2))
       end if
       
       return
       end


       function efficiency(risc)
c      Calculates the accretion efficiency, s.t L_bol=eta*Mdot*c^2
c      Uses GR case where eta = 1 - sqrt(1 - 2/(3*risco))
c      Taken from: The Physics and Evolution of Active Galactic Nuclei,
c      H. Netzer, 2013, p.38
       
       implicit none
       double precision risc, efficiency

       efficiency = 1.0 - sqrt(1.0 - 2.0/(3.0*risc))
       return
       end
       

       function calc_rsg(M, mdot)
c      Calculates the self gravity radius following Laor & Netzer 1989
c
c      Note! Assumes \alpha=0.1

       implicit none
       double precision M, mdot, calc_rsg
       double precision m9, alpha

       alpha = 0.1
       m9 = M/1.0d9
       calc_rsg = 2150*m9**(-2.0/9.0)*mdot**(4.0/9.0)*alpha**(2.0/9.0)

       return
       end


       function NTpars(r, astar, risc)
c      Calculates the Novikov-Thorne parameters at a given radius

       implicit none
       double precision r, astar, risc
       double precision pi, y, yisc, y1, y2, y3
       double precision B, C1, C2, C, NTpars
       double precision C2_1, C2_2, C2_3

       pi = 4.0*atan(1.0)

       y = sqrt(r)
       yisc = sqrt(risc)
       y1 = 2.0*cos((acos(astar) - pi)/3.0)
       y2 = 2.0*cos((acos(astar) + pi)/3.0)
       y3 = -2.0*cos(acos(astar)/3.0)

       B = 1.0 - (3.0/r) + ((2.0*astar)/(r**(3.0/2.0)))

       C1 = 1 - (yisc/y) - ((3.0*astar)/(2.0*y)) * log(y/yisc)

       C2_1 = (3.0*(y1-astar)**2.0)/(y*y1 * (y1-y2) * (y1-y3))
       C2_1 = C2_1 * log((y-y1)/(yisc-y1))
       
       C2_2 = (3.0*(y2-astar)**2.0)/(y*y2 * (y2-y1) * (y2-y3))
       C2_2 = C2_2 * log((y-y2)/(yisc-y2))

       C2_3 = (3.0*(y3-astar)**2.0)/(y*y3 * (y3-y1) * (y3-y2))
       C2_3 = C2_3 * log((y-y3)/(yisc-y3))

       C2 = C2_1 + C2_2 + C2_3
       C = C1 - C2
       NTpars = C/B

       return
       end


       function nttemp(r, M, mdot, Mdot_edd, astar, risc)
c      Function to calculate Novikov-Thorne temperature at radius r
      
       implicit none

       double precision r, M, mdot, Mdot_edd, astar, risc
       double precision pi, G, sigma_sb, c, Rg, Rt
       double precision NTpars
       double precision nttemp

       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33 !cm^3 s^-1 Msol^-1
       sigma_sb = 5.67d-5 !erg cm^-2 s^-1 K^-4
       c = 3.0d10 !cm/s

       Rg = (G*M)/c**2 !cm
       Rt = NTpars(r, astar, risc)

       nttemp = (3.0*G*M*mdot*Mdot_edd)/(8*pi*sigma_sb*(r*Rg)**3) !K^4
       nttemp = nttemp*Rt
       nttemp = nttemp**(1.0/4.0) !K

       return
       end
       
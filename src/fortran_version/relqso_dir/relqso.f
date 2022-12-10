c      Updated version of QSOSED
c      Ref. Hagen and Done (in prep.)
c
c      The model consists of two regions: a warm
c      Comptonization region where the disc has not thermalised properly,
c      and a spherical inner corona modelled as a hot Comptonization
c      region. All the energetics are calculated in the rest frame of the
c      black hole, before the relativistic transfer functions are applied
c      The Comptonization is modelled using NTHCOMP (Zdziarski, Johnson &
c      Magdziarz 1996; Zycki, Done & Smith 1999)
c
c      The disc truncation radius (Rhot) is set by the requirement that
c      the dissipated luminosity in the corona Lx_diss=0.02Ledd
c      This then also give Gamma_hot from Eqn. 6 in Kubota and Done (2018)
c
c      Additionaly, this subroutine includes the option for relativistic
c      effects to be taken into account, using KYCONV (Dovciak et al. 2004)
      
       subroutine relqso(ear,ne,param,ifl,photar,photer)
       implicit none

       integer npars
       parameter(npars=7)
      
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

       integer i, n             !Iteration indeces
       double precision zfac

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
c      param(6):  redshift
c      param(7):  rel, flag, 1=relativistic, 0=non-relativistic

c      checking if parameters have changed
       parchange=.false.
       do i=1,npars,1
          if (param(i).ne.oldpar(i)) parchange=.true.
       end do
       
c      Checking if we need to change energy grid
c      The default extends from 1e-4 to 1e3 keV
       zfac = (1.0 + param(6))
       if (ear(0).eq.0.0) then
          newemin=ear(1) - (ear(2)-ear(1))/10.0
       else
          newemin=min(1.0e-4, ear(0)*zfac)
       end if
       newemax=max(1.0e3, ear(ne)*zfac)

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

c      Call model if parameters or energy grid have changed
       if (parchange.or.echange) then
          
          call calc_relqso(enew, Nnew, param, ifl, ph)

          !Redshift correct energy bins
          do n=1, Nnew, 1
             enew(n) = enew(n)/zfac
             ph(n) = ph(n)/zfac
          end do

          do i=1, npars, 1
             oldpar(i)=param(i)
          end do
       end if

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


       subroutine calc_relqso(es, nn, param, ifl, ph)

c      Function to calculate the SED. Calculates the spectral
c      contribution from each region (warm comp, hot comp)
c      Temperature found through the Novikov-Thorne equiations
c       
c      The relativistic corrections are applied at each annulus through
c      kyconv

       implicit none
       integer nn, ifl, rel, show_comp
       real es(0:nn), ph(nn), param(*)
       real ph_ann(nn), phw_ann(nn), phh_ann(nn), phh_shape(nn) !disc, warm, hot
       double precision M, mdot
       double precision rout, rw, rh, risco
       double precision astar, hmax
       double precision kTh, kTw
       double precision gammah, gammaw
       double precision cos_inc, dist

       double precision pi, G, c, h, sigma_sb,  mp, kB
       double precision eta, Ledd, Mdot_edd, Rg, Ldiss, Lseed
       double precision Lwarm, Ldisc
       double precision dr_dex, nrd, nrw, nrh, nrseed, nrall, dex_itr
       double precision dlog_rd, dlog_rw, dlog_rh, dlog_rseed, dlog_rall
       double precision rmid, dr
       double precision kkev, kevhz
       double precision dflux, en, dnphot, tr
       double precision calc_risco, calc_rsg, calc_fcol
       double precision efficiency, nttemp, reptemp
       double precision lphseed, fcov, theta0
       double precision lseed_ann, lseed_cann, ldiss_ann, ldiss_itr
       double precision tseed_h, ysb, Ldiss_act

       real kpar(12), ntwpar(5), nthpar(5)
       real kpherr(nn), ntwpherr(nn), nthpherr(nn)

       double precision normw_ann, normh_ann
       double precision ntw_out, nth_out

       logical return_warm, return_hot, return_disc
       logical hotlim, userel
       integer i, j, n

       double precision ltest

c-----------------------------------------------------------------------
c      Section for setting up before main caluclation
c-----------------------------------------------------------------------
c      Constants
       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33   !Grav const. cm^-3 s^-1 Msol^-1
       c = 3.0d10              !Speed of light. cm/s
       h = 6.62617d-27         !Plank const. erg s
       sigma_sb = 5.670374d-5  !Stefan-Boltzmann const. erg s^-1 cm^-2 K^-4
       mp = 1.67d-24           !Proton mass. g
       kB = 1.38d-16           !Boltzmann const. erg/K

c      Unit conversion constants
       kkev = 1.16048d7          !K/keV
       kevhz = 2.417965d17      !Hz/keV

       dr_dex = 50              !Radial resolution, Nbins per decade (ish)
      

c      Read parameters
       M = dble(param(1))
       dist = dble(param(2)) !Mpc
       mdot = dble(10**(param(3)))
       astar = dble(param(4))
       cos_inc = dble(param(5))
       rel = param(7)

       if (rel.eq.0) then
          userel=.false.
       else
          userel=.true.
       end if
       
       
c      Getting accretion attributes
       risco = calc_risco(astar)
       eta = efficiency(risco)
       rout = calc_rsg(M, mdot) !Outer radius always self-gravity
       
       Ledd = 1.39d38 * M       !erg/s
       Ldiss = 0.02 * Ledd !erg/s
       Mdot_edd = Ledd/(eta * c**2) !g/s
       Rg = (G*M)/c**2          !cm
       
       dist = dist*1.0d6        !pc
       dist = dist*3.086d18     !cm

c      Setting and calculating Compton params
       kTh = dble(100.0)
       kTw = dble(0.2)
       gammaw = 2.5
       

c      Calculating Rhot, based off requirement Ldiss = 0.02Ledd
c      Uses a more refined grid to avoind innacuracies in inner region
c      NOTE: Requirement set in rest frame!!!!
       nrall = (log10(rout) - log10(risco))*1000 !Nr bins over entire flow
       nrall = ceiling(nrall)    !round up to nearest integer
       dlog_rall = (log10(rout) - log10(risco))/nrall !actual bin size
       ldiss_ann = dble(0.0)    !Initialising
       ldiss_itr = dble(0.0)
       hotlim = .true.

       ditr: do i=1, int(nrall), 1
          rmid =10**(log10(risco)+float(i-1)*dlog_rall+dlog_rall/2.0)
          dr = 10**(log10(rmid)+dlog_rall/2.0)
          dr = dr - 10**(log10(rmid)-dlog_rall/2.0)
 
          if ((rmid+dr/2.0 + dr).gt.rout) then !Makes big bin at edge
             dr = rout - 10**(log10(rmid)-dlog_rall/2.0)
          end if

          tr = nttemp(rmid, M, mdot, Mdot_edd, astar, risco) !Temp. 
          ldiss_ann = ldiss_ann+4.0*pi*rmid*dr*Rg**2.0*sigma_sb*tr**4.0 !erg/s
          ldiss_ann = sigma_sb * tr**4.0
          ldiss_ann = ldiss_ann * 4.0*pi*rmid*dr*Rg**2.0

          ldiss_itr = ldiss_itr + ldiss_ann

          rh = rmid
          if (ldiss_itr.gt.Ldiss) then
             write(*,*) ldiss_itr/Ldiss
             hotlim = .false.
             exit ditr               !break condition
          end if
       end do ditr

       rw = 2*rh
       hmax = dble(min(10.0, rh)) !Coronal scale height


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
       ntwpar(4) = 0.0 !BB assumed from each annulus
       ntwpar(5) = 0.0 !Redshift - dealt with within this script

c      And finally hot nthomp params
       nthpar(2) = kTh !Note gammah set later after calculating
       nthpar(4) = 0.0
       nthpar(5) = 0.0


       return_disc=.true.
       return_warm=.true.
       return_hot=.true.
       
c      Checking we actually reached Ldiss = 0.02Ledd
       if (hotlim) then
          call xwrite('WARNING!!! Ldiss never reaches 0.02Ledd', 5)
          call xwrite('=> No upper limit for r_hot', 5)
          call xwrite('Try increasing mass accretion rate!', 5)
          return_warm=.false.
          return_disc=.false.
          rh = rout
       end if

c      Checking rw !> rout
       if (rw.gt.rout) then
          call xwrite('WARNING!!! rw > rout ---- Setting rw = rout', 10)
          return_disc=.false.
       end if



c-----------------------------------------------------------------------
c      Section for calculating the disc region
c-----------------------------------------------------------------------
       if (return_disc) then
          nrd = (log10(rout) - log10(rw)) * dr_dex !nr of bins in region
          nrd = ceiling(nrd) !round up to nearest integer value
          dlog_rd = (log10(rout) - log10(rw))/nrd
          Ldisc = dble(0.0) !Re-setting disc lum
          do i=1, int(nrd), 1
             rmid = 10**(log10(rw)+float(i-1)*dlog_rd+dlog_rd/2.0)
             dr = 10**(log10(rmid)+dlog_rd/2.0)
             dr = dr - 10**(log10(rmid)-dlog_rd/2.0)

             if ((rmid+dr/2.0 + dr).gt.rout) then !Makes big bin at edge
                dr = rout - 10**(log10(rmid)-dlog_rd/2.0)
             end if

             tr = reptemp(rmid,M,mdot,Mdot_edd,astar,risco,hmax,Ldiss) !Temp. K
             
c            Calculating BB emission over each energy
             do n=1, nn, 1
                en = dble(log10(es(n-1)) + log10(es(n)))
                en = en/2.0 !Evaluate in mid bin
                en = 10**en
                
                if (en.lt.30.0*tr) then
                   dflux = (pi*2.0 * h * (en*kevhz)**3.0)/(c**2.0)
                   dflux = dflux * 1/(exp((h*en*kevhz)/(kB*tr)) - 1.0)
                   dflux = dflux * 4.0 *  Rg**2.0 !Lum in erg/s/Hz
                   
                   dflux = dflux/(4.0*pi*dist**2.0) !erg/s/cm^2/Hz
                   
                   ph_ann(n) = dflux/(h*kevhz*en) !photons/s/cm^2/Hz
                   ph_ann(n) = ph_ann(n)*kevhz*(es(n)-es(n-1)) !photons/s/cm^2                
                else
                   ph_ann(n) = 0.0d0
                end if         
             end do

c            Applying kyconv
             if (((rmid + dr/2.0).lt.1.0d3).and.(userel)) then
                kpar(3) = rmid - dr/2.0
                kpar(5) = rmid + dr/2.0
                kpar(8) = rmid
                call kyconv(es, nn, kpar, ifl, ph_ann, kpherr)
             else
                ph_ann = ph_ann * pi*rmid*dr*(cos_inc/0.5)
             end if

             Ldisc = Ldisc + 4.0*pi*rmid*dr*Rg**2*sigma_sb*tr**4.0 !erg/s

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

c-----------------------------------------------------------------------
c      Section for calculating the warm Compton region
c-----------------------------------------------------------------------
       if (return_warm) then
          nrw = (log10(rw) - log10(rh)) * dr_dex !nr bins in region
          nrw = ceiling(nrw)    !round up to nearest integer
          dlog_rw = (log10(rw) - log10(rh))/nrw !actual bin size
          Lwarm = dble(0.0)
          do i=1, int(nrw), 1
             rmid = 10**(log10(rh)+float(i-1)*dlog_rw + dlog_rw/2.0)
             dr = 10**(log10(rmid) + dlog_rw/2.0)
             dr = dr - 10**(log10(rmid) - dlog_rw/2.0)

             if ((rmid+dr/2.0+dr).gt.rw) then !Ensuring bin within region
                dr = rw - 10**(log10(rmid) - dlog_rw/2.0)
             end if

             tr = reptemp(rmid,M,mdot,Mdot_edd,astar,risco,hmax,Ldiss) !Temp. K
             ntwpar(3) = tr/kkev !setting seed temp to disc annulus temp (keV)

c            Calling nthcomp
             call donthcomp(es, nn, ntwpar, ifl, phw_ann, ntwpherr) !photons/s/cm^s/

             normw_ann = sigma_sb*tr**4.0
             normw_ann = normw_ann * 4.0*Rg**2.0 !erg/s/, emission from annulus
             normw_ann = normw_ann/(4.0*pi*dist**2.0) !erg/s/cm^2/
c            Again leaving off pi*r*dr * cosi/0.5 as this done in kyconv!!!
             
c            Finding total flux output from nthcomp
             ntw_out = 0.0
             do n=1, nn, 1
                ntw_out = ntw_out + phw_ann(n)*es(n)*kevhz*h !erg/s/cm^2
             end do

c            Now applying normalisation
             do n=1, nn, 1
                if (ntw_out.eq.0) then
                   phw_ann(n) = 0.0
                else
                   phw_ann(n) = phw_ann(n) * (normw_ann/ntw_out) !photons/s/cm^2
                end if
             end do

c            Applying kyconv
             if (((rmid+dr/2.0).lt.1.0d3).and.(userel)) then
                kpar(3) = rmid - dr/2.0
                kpar(5) = rmid + dr/2.0
                kpar(8) = rmid
                call kyconv(es, nn, kpar, ifl, phw_ann, kpherr)
             else
                phw_ann = phw_ann * pi*rmid*dr*(cos_inc/0.5)
             end if

             Lwarm = Lwarm + 4.0*pi*rmid*dr*Rg**2*sigma_sb*tr**4.0 !erg/s
            

c            Adding to total output array
             do n=1, nn, 1
                if (i.eq.1) then
                   if (return_disc) then
                      ph(n) = ph(n) + phw_ann(n)
                   else
                      ph(n) = phw_ann(n) !start from scratch
                   end if
                else
                   ph(n) = ph(n) + phw_ann(n)
                end if
             end do     
          end do
       end if



c-----------------------------------------------------------------------
c      Section for calculating hot compton region (ie corona)
c-----------------------------------------------------------------------
       if (return_hot) then
c         First finding the total seed photon luminosity/flux
c         Integrated flux from entire disc seen by corona
          nrseed = (log10(rout) - log10(rh)) * dr_dex
          nrseed = ceiling(nrseed) !round up to nearest integer
          dlog_rseed = (log10(rout) - log10(rh))/nrseed

          lphseed = 0.0
          do i=1, int(nrseed), 1
             rmid = 10**(log10(rh)+float(i-1)*dlog_rseed+dlog_rseed/2.0)
             dr = 10**(log10(rmid) + dlog_rseed/2.0)
             dr = dr - 10**(log10(rmid) - dlog_rseed/2.0)

             if ((rmid+dr/2.0+dr).gt.rout) then !Ensuring bin within region
                dr = rout - 10**(log10(rmid) - dlog_rseed/2.0)
             end if

             if (hmax.le.rmid) then
                theta0 = asin(hmax/rmid)
                fcov = theta0 - 0.5*sin(2.0*theta0) !corona covering fraction seen from rmid
             else
                fcov = 0.0
             end if

             tr = nttemp(rmid, M, mdot, Mdot_edd, astar, risco) !Temp. K

             lseed_ann = sigma_sb * tr**4.0 !erg/s/cm^2
             lseed_ann = lseed_ann * 4.0*pi*rmid*dr*Rg**2.0 * (fcov/pi) !erg/s
             lseed_ann = lseed_ann/(4*pi * dist**2.0) !erg/s/cm^2

             lphseed = lphseed + lseed_ann !total seed photon flux erg/s/cm^2
          end do

          Lseed = lphseed * 4*pi * dist**2.0 !erg/s
          gammah = (7.0/3.0) * (Ldiss/Lseed)**(-0.1)
          nthpar(1) = gammah

c         Now finding seed photon temperature
c         Assumed to be ~inner disc temp
          tseed_h = nttemp(rh, M, mdot, Mdot_edd, astar, risco) !K
          ysb = (gammaw*(4.0/9.0))**(-4.5) !Compton y-param for warm region
          tseed_h = (tseed_h/kkev) * exp(ysb) !in keV
         
          
          !calling nthcomp now, since intrinsic shape does not change!!!
          nthpar(3) = tseed_h
          call donthcomp(es, nn, nthpar, ifl, phh_shape, nthpherr) 
          !total flux output from nthcom
          nth_out = 0.0
          do n=1, nn, 1
             nth_out = nth_out + phh_shape(n)*es(n)*kevhz*h !erg/s/cm^2
          end do

c         Now calculating emission from each coronal annulus
          nrh = (log10(rh) - log10(risco)) * dr_dex
          nrh = ceiling(nrh)    !round up to integer
          dlog_rh = (log10(rh) - log10(risco))/nrh !actual spacing

          ltest = dble(0.0)
          do i=1, int(nrh), 1
             rmid = 10**(log10(risco)+float(i-1)*dlog_rh+dlog_rh/2.0)
             dr = 10**(log10(rmid) + dlog_rh/2.0)
             dr = dr - 10**(log10(rmid) - dlog_rh/2.0)

             if ((rmid+dr/2.0+dr).gt.rh) then !Ensuring bin within region
                dr = rh - 10**(log10(rmid) - dlog_rh/2.0)
             end if

             tr = nttemp(rmid, M, mdot, Mdot_edd, astar, risco) !K
             ldiss_ann = sigma_sb * tr**4.0
             ldiss_ann = ldiss_ann * 4.0 * Rg**2.0 !erg/s

             ltest = ltest + ldiss_ann * pi*rmid*dr
             ldiss_ann = ldiss_ann/(4.0*pi * dist**2.0) !erg/s/cm^2
             !Assuming equall amounts of emission from each coronal annulus
             !Need to renorm here for kyconv
             lseed_cann = lphseed/(pi*rmid*dr * nrh)

             normh_ann = ldiss_ann + lseed_cann !erg/s/cm^2

             !Now applying normalisation
             do n=1, nn, 1
                if (nth_out.eq.0) then
                   phh_ann(n) = 0.0
                else
                   phh_ann(n) = phh_shape(n) * (normh_ann/nth_out) !photons/s/cm^2
                end if
             end do

c            Applying kyconv
             if (((rmid+dr/2.0).lt.1.0d3).and.(userel)) then
                kpar(3) = rmid - dr/2.0
                kpar(5) = rmid + dr/2.0
                kpar(8) = rmid
                call kyconv(es, nn, kpar, ifl, phh_ann, kpherr)
                phh_ann = phh_ann/(2.0*cos_inc) !not disc so rmoving cos(i) dependence
             else
                phh_ann = phh_ann * pi*rmid*dr
             end if

c            Adding to total output array
             do n=1, nn, 1
                if (i.eq.1) then
                   if ((return_warm).or.(return_disc)) then
                      ph(n) = ph(n) + phh_ann(n)
                   else
                      ph(n) = phh_ann(n) !if no warm comp start from scratch
                   end if
                else
                   ph(n) = ph(n) + phh_ann(n)
                end if
             end do
          end do
       end if

c     Writing outpur parameters IF wanted
       write (*,*) 'Disc Pars'
       write (*,*) '-----------------------------------------'
       write (*,*) 'r_hot = ',rh,'rw = ',rw,'r_out = ',rout
       write (*,*) '-----------------------------------------'
       write (*,*)
       write (*,*) 'Luminosities'
       write (*,*) '-----------------------------------------'
       write (*,*) 'Ldiss/Ledd = ',ltest/Ledd
       write (*,*) 'Lhot/Ledd = ',(Ldiss+Lseed)/Ledd
       write (*,*) 'Lhot=',Ldiss+Lseed,'erg/s'
       write (*,*) 'Lwarm=',Lwarm,'erg/s'
       write (*,*) 'Ldisc=',Ldisc,'erg/s'
       write (*,*) 'Ltot=',Ldiss+Lseed+Lwarm+Ldisc,'erg/s'
       write (*,*) '-----------------------------------------'
       write (*,*)
       write (*,*) 'Spec Pars'
       write (*,*) '-----------------------------------------'
       write (*,*) 'gammah = ',gammah,' , gammaw = ',gammaw
       write (*,*) 'kT_hot_seed = ',tseed_h,'keV , y_param = ',ysb
       write (*,*)
      
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
       sigma_sb = 5.670367d-5 !erg cm^-2 s^-1 K^-4
       c = 3.0d10 !cm/s

       Rg = (G*M)/c**2.0 !cm
       Rt = NTpars(r, astar, risc)

       nttemp = (3.0*G*M*mdot*Mdot_edd)/(8*pi*sigma_sb*(r*Rg)**3) !K^4
       nttemp = nttemp*Rt
       nttemp = nttemp**(1.0/4.0) !K

       return
       end


       function reptemp(r, M, mdot, Mdot_edd, astar, risc, hmax, Lx)
c      Disc temperature inclusing re-processing

       implicit none

       double precision r, M, mdot, Mdot_edd, astar, risc, hmax, Lx
       double precision pi, G, sigma_sb, c, Rg, Rt
       double precision reptemp, nttemp, ftemp
       double precision frep, Rphys, Hphys

       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33 !cm^3 s^-1 Msol^-1
       sigma_sb = 5.670367d-5 !erg cm^-2 s^-1 K^-4
       c = 3.0d10               !cm/s

       Rg = (G*M)/c**2.0  !cm
       Rphys = r*Rg !cm
       Hphys = hmax*Rg !cm
       
       frep = (0.5*Lx)/(4.0*pi * (Rphys**2.0 + Hphys**2.0)) !erg/s/cm^2
       frep = frep * Hphys/sqrt(Rphys**2.0 + Hphys**2.0)
       frep = frep * (1.0-0.3)    !Disc albedo assumed as 0.3

       ftemp = frep/sigma_sb

       reptemp = nttemp(r, M, mdot, Mdot_edd, astar, risc)
       reptemp = reptemp**4.0
       reptemp = reptemp + ftemp
       reptemp = reptemp**(1.0/4.0) !K

       return
       end


        

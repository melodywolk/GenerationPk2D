import pyfits
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import fnmatch
import os
import cosmolopy.distance as cp
import cosmolopy.parameters as p
import pt
import scipy
from scipy.integrate import simps
from scipy import interpolate


def projection(Nofz,z, biasz, Nboxsize ,Volume, Nli, Poisson, Window):
    print "BIAS"
    print biasz
    # Multipole parameters : log grid from lmin to lmax with Nl points
    Nl = Nli
    lmin = 20.
    lmax = 8000.
    lnl = np.linspace(np.log(lmin),np.log(lmax),Nl)
    l=np.exp(lnl)
    # Numbers of bins :
    Nbins = len(Nofz[:,0])
    ztab = np.array(Nofz[:,0])
    dndz = np.array(Nofz[:,1])
    # Main cosmological parameters
    h = 0.70
    omb = 0.045
    omega_m = 0.27
    omega_l = 1-omega_m
    n_fid=0.96
    sig8 =0.80
    omnu = 0
    N_nu = 0
    omc = omega_m - omb - omnu
    cosmo={'omega_M_0' : omega_m, 'omega_lambda_0' : omega_l, 'h' : h, 'w' : -1, 'omega_n_0' : omnu, 'N_nu' : N_nu}
    print cosmo
    # dchi/dz = c/H, in Mpc
    cosmo = cp.set_omega_k_0(cosmo)
    dchidz = np.zeros((Nbins))
    chi = np.zeros((Nbins))
    for i in range(0,Nbins):
        dchidz[i] = cp.comoving_integrand(ztab[i], **cosmo)
        chi[i] = cp.comoving_distance(ztab[i], **cosmo)
        #print "Comoving distance at z=%lf is %lf Mpc and dchidz=%lf" % (ztab[i], chi[i], dchidz[i])
    
    # 3D power spectrum 
    # Call CAMB for the power spectrum, uses HALOFIT
    lnkmax = np.log(lmax) - np.log(np.min(chi)) # maximal 3d lnk needed, k 3d  in 1/Mpc 
    c = pt.Camb(hubble=h*100.,ombh2=omb*h**2,transfer_redshift = 0,transfer_kmax = np.exp(lnkmax)/h,do_nonlinear=0,omega_lambda=omega_l,omch2=omc*h**2) 
    c.run()
    R=8 # Mpc
    sigz2 = pt.sigma2fromPk(c,R)
    print 'sigma_8'
    print np.sqrt(sigz2)
    ampl = c.cp.scalar_amp[0]
    f = float((sig8**2)/sigz2)
    c = pt.Camb(hubble=h*100.,ombh2=omb*h**2,transfer_redshift=[z],transfer_kmax=np.exp(lnkmax)/h,do_nonlinear=1,omega_lambda=omega_l,omch2=omc*h**2,scalar_amp = [ampl*f])
    c.run()
    power = np.loadtxt(open('tmpCosmoPyCamb_matterpower.dat')) #Pk has units fof (Mpc/h)^3
    k = np.array(power[:,0])
    Pk = np.array(power[:,1]) 

    np.savetxt('verif.txt', np.transpose([k,Pk])) 

    lnk = np.log(k) + np.log(h) #in 1/Mpc
    k = np.exp(lnk) #1/Mpc
    lnpk = np.log(Pk) - 3.0*np.log(h) # lnpk better than Pk for interpolation. Now in Mpc^3
    Nk = len(lnk)
    print 'kmax in h/Mpc %lf' % (np.exp(lnkmax)/h)
    # Interpolates the 3d power spectrum P(k,z) at the requested arguments (l/chi(z), z)
    Pk_zl = np.zeros((Nl)).reshape((Nl))
    Pk_zl = np.exp(np.interp(lnl-np.log(np.mean(chi)), lnk, lnpk))
    # Now we calculate the Cl for each cell
    Cl = np.zeros((Nl)).reshape((Nl))
    chi2dchidz = chi*chi*dchidz
    pz=dndz/simps(dndz,ztab)
    for i in range(0, Nl):
	Cl[i] = simps(pz*pz/chi2dchidz*Pk_zl[i], ztab)*biasz*biasz

    if(Window==1):
     # Pixel window functions for each auto and cross spectra if volume and boxsize are provided :
        kl = l*np.sqrt(Volume)/(2.*Nboxsize)
        xp=np.arange(0,np.pi, 0.001)
        for i in range(0, Nl):
            Cl[i] *= (1./np.pi)*simps(np.sinc((kl[i]*np.cos(xp))/(np.pi))**2*np.sinc((kl[i]*np.sin(xp))/(np.pi))**2, xp)

    if(Poisson==1):
     #  Shot noise if volumes are provided requested :
        for i in range(0,Nl-1):
            Cl[i] += 1./ ( np.sum(Nofz[:,1]) /Volume)

    return l[:-1], Cl[:-1]
    

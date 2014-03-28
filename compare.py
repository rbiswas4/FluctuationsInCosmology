#!/usr/bin/env python
#Comparison of power spectrum 
#	from the CAMB power spectrum and CAMB transfer function output
#
#Requirements :
#	(a) Make sure that Suman's code directory is available
#	(b) Make sure that interfacecosmology is available
#Data Files:
#    T	
#Comparison with Suman's code 
#Ecample with m000 to check that things are working
import numpy as np
import matplotlib.pyplot as plt
import camb_utils.cambio as cio
import psutils as psu
import fsigma as fs
import MF as MF
from interfaces import FCPL
import massfunctions as mf
import plotutils as pu 
import sys

# INPUT FILES 
dirn  = "interfacecosmology/"
tkfile  = dirn +"example_data/M000/m000n0_transfer_fin_out.dat"
pkfile = dirn + "example_data/M000/m000n0_matterpower_fin.dat"
#psrunfile = "/home/rbiswas/doc/CAMB_outs/m000r03/m000r03_matterpower_fin.dat"

M000 = FCPL(H0 = 71. ,Om0  = 0.2647, Ob0 = 0., ns = 0.963, As = 2.16e-9)
M000s = FCPL(H0 = 71. ,Om0  = 0.2647, Ob0 = 0., ns = 0.963, sigma8 = 0.8)
M000p = FCPL(H0 = 71. ,Om0  = 0.3647, Ob0 = 0.0, ns = 0.963, As = 2.16e-9)

Tks = cio.loadtransfers(filename = tkfile)
k = Tks[:,0]
Tk = Tks[:,-1]
sumanrhoc    = 2.77536627e11
sumannorm    = fs.Pk_norm(sigma8 = 0.8, ns = 0.963, k = k , Tk=Tk, hubble= 0.71)

	#working power spectrum example
psfrompk = psu.powerspectrum(koverh = None, asciifile =pkfile)
	#working transfer function example
psfromtk = psu.powerspectrum(koverh = None, pstype ="cb", asciifile =tkfile, cosmo = M000)
psfromtkm = psu.powerspectrum(koverh = None, pstype ="matter", asciifile =tkfile, cosmo = M000)
psfromtkms = psu.powerspectrum(koverh = None, pstype ="matter", asciifile =tkfile, cosmo = M000s)
psfromtkcbs = psu.powerspectrum(koverh = None, pstype ="cb", asciifile =tkfile, cosmo = M000s)
psfromtkcbscb = psu.powerspectrum(koverh = None, pstype ="cb",sigma8type= "cb", asciifile =tkfile, cosmo = M000s)
#psfromtkcbscb = psu.powerspectrum(koverh = None, pstype ="cb", asciifile =tkfile, cosmo = M000s)
#print  len(psfromtk), "now"
pscompfig, pscomp_ax0, pscomp_ax1 = pu.settwopanel(setdifflimits = [0.99,1.01])
#pscompfig, pscomp_ax0, pscomp_ax1 = pu.settwopanel(setdifflimits = None)
pscomp_ax0.loglog(psfromtk[0], psfromtk[1], label = "TKcb As=2.16e-9")
pscomp_ax0.loglog(psfromtkm[0], psfromtkm[1], label = "TKm, As = 2.16e-9")
pscomp_ax0.loglog(psfromtkms[0], psfromtkms[1], label = "TKm, sigma8 =0.8")
pscomp_ax0.loglog(psfromtkcbs[0], psfromtkcbs[1], label = "TKcb, sigma8 =0.8")
pscomp_ax0.loglog(psfromtkcbscb[0], psfromtkcbscb[1], label = "TKcb, sigma8 cb =0.8")
pscomp_ax0.loglog(psfrompk[0], psfrompk[1], label = "PK")
pscomp_ax0.legend(loc="best")

pscomp_ax1.plot(psfromtk[0], psfromtk[1]/psfrompk[1],label = "Tkcb, As = 2.16e-9")
#pscomp_ax1.plot(psfromtk[0], psfromtkm[1]/psfrompk[1],label = "Tkm, As = 2.16e-9")
pscomp_ax1.plot(psfromtk[0], psfromtkcbs[1]/psfrompk[1],label = "Tkcb, sigma8 = 0.8")
pscomp_ax1.plot(psfromtk[0], psfromtkcbscb[1]/psfrompk[1],label = "Tkcb, sigma8 cb = 0.8")
pscomp_ax1.legend(loc = "best")
#pscomp_ax0.set_xscale('log')
#print psfrompk[1]/psfromtk[1]



#Masses in units of solar masses
Masses = np.logspace(8,16,100)
Rad = np.linspace(.1,20,100)
MassesinhoverM  = Masses* M000.h
sumansigmar = np.zeros(len(Rad))
sumandlsidlm = np.zeros(len(Masses))
sumandndlnM = np.zeros(len(Masses))
sumansigmam = np.zeros(len(Masses))
sumanfsigma = np.zeros(len(Masses))
for i,rad in enumerate(Rad) :
	s = fs.sigmar ( Om = M000.Om0, ns = 0.963, R = rad, N=sumannorm , k = k , Tk = Tk , hubble = 0.71 )
	#print rad , s
	sumansigmar[i] = s
	#print sumansigmar
for i, mass in enumerate(MassesinhoverM) :
	sumansigmaM  = fs.sigmam ( Om = M000.Om0, ns = 0.963, M= mass, N=sumannorm , k = k , Tk = Tk , hubble = 0.71 ) 
	sumanlogsigmaM = fs.logsigm ( Om = M000.Om0 , ns = 0.963, M= mass, N = sumannorm, k = k, Tk = Tk, sigmam = sumansigmaM, hubble = 0.71 )
	sumanfsigma[i] = MF.MF_fit(sumansigmaM , z = 0.)
	sumandndlnM [i] = sumanfsigma [i]*(sumanrhoc*M000.Om0*sumanlogsigmaM)/mass
	#sumandndlnM [i] = sumanfsigma [i]*(sumanlogsigmaM)
	sumansigmam [ i] = sumansigmaM 
	sumandlsidlm [ i] = sumanlogsigmaM 

massfnpkM000 = psu.dndlnM(M= Masses 
	, ps = psfrompk
	, cosmo = M000)
#massfnpkM000p = psu.dndlnM(M= Masses , ps = psfrompk, cosmo = M000p)

	#sigma as a function of distance
sigmarfig, sigmar_ax0 , sigmar_ax1 = pu.settwopanel()
sigmar_ax0.plot(Rad, sumansigmar, "o",label = "Suman TF")
sigmar_ax0.plot(Rad, psu.sigma(R = Rad, ps = psfrompk, cosmo = M000),"+", label = "PS")
sigmar_ax0.set_ylabel(r'$\sigma (R) $')
sigmar_ax0.legend(loc= "best", numpoints =1)
sigmar_ax1.plot(Rad, sumansigmar/ psu.sigma(R = Rad, ps = psfrompk, cosmo = M000),"+", label = "sumanTF/PS")
sigmar_ax1.set_xlabel(r'$R \/(h^{-1}Mpc) $')

	#sigma as a function of halo mass
sigmaMfig, sigmaM_ax0 , sigmaM_ax1 = pu.settwopanel()
sigmaM_ax0.plot(MassesinhoverM , sumansigmam, "o",label = 'Suman TF')
sigmaM_ax0.plot(MassesinhoverM, psu.sigmaM(M= Masses, ps = psfrompk, cosmo = M000),"+", label = 'PS')
sigmaM_ax0.set_xscale('log')
sigmaM_ax0.set_ylabel(r'$\sigma(M)$')
sigmaM_ax0.legend(loc= "best", numpoints = 1)
sigmaM_ax1.plot(MassesinhoverM ,sumansigmam/psu.sigmaM(M= Masses, ps = psfrompk, cosmo = M000))
sigmaM_ax1.grid(True)
sigmaM_ax1.set_xscale('log')
sigmaM_ax1.set_xlabel(r'$M (h^{-1}M_\odot) $')

	#fsigma as a function of halo mass
dfsigfig , dfsig_ax0, dfsig_ax1 = pu.settwopanel()
dfsig_ax0.plot(MassesinhoverM , sumanfsigma , "g", label = "Suman TF") 
dfsig_ax0.plot(MassesinhoverM , mf.__fsigmaBhattacharya(psu.sigmaM(M = Masses, ps = psfrompk, cosmo = M000)) , "rs", label = "PS") 
dfsig_ax0.set_ylabel(r'$f_{\sigma}(M)$')
dfsig_ax0.legend(loc  = "best", numpoints = 1)
dfsig_ax0.set_xscale('log')
dfsig_ax1.plot(MassesinhoverM , sumanfsigma/ mf.__fsigmaBhattacharya(psu.sigmaM(M = Masses, ps = psfrompk, cosmo = M000)) , "rs", label = "PS") 
dfsig_ax1.set_xscale('log')

	#dlnsigmainv dlnM as a function of halo mass
dlsidlMfig , dlsidlM_ax0, dlsidlM_ax1 = pu.settwopanel()
dlsidlM_ax0.loglog( MassesinhoverM , -sumandlsidlm, "go", label = "TF") 
dlsidlM_ax0.loglog( MassesinhoverM , psu.dlnsigmadlnM(M= Masses, ps = psfrompk, cosmo = M000) , "r+" , label = "PS") 
dlsidlM_ax0.set_ylabel(r'$dlnsigma/dlnM$')
dlsidlM_ax0.legend(loc = "best", numpoints = 1)
dlsidlM_ax1.plot( MassesinhoverM , sumandlsidlm / psu.dlnsigmadlnM (M = Masses,  ps = psfrompk, cosmo = M000))
dlsidlM_ax1.grid(True)
dlsidlM_ax1.set_xscale("log")

	#mass function as a function of halo mass

#dndlnmfig, dndlnm_ax0 , dndlnm_ax1 = pu.settwopanel(setdifflimits = None)
dndlnmfig, dndlnm_ax0 , dndlnm_ax1 = pu.settwopanel()
dndlnm_ax0.loglog(MassesinhoverM ,  -sumandndlnM , label = "Suman")
dndlnm_ax0.loglog(MassesinhoverM , psu.dndlnM(M = Masses, ps = psfrompk, cosmo = M000)/M000.h**3.0, label = "PS")
dndlnm_ax0.set_ylabel(r'$dn/d\ln{(M)} (h^3 Mpc^{-3} $')
dndlnm_ax0.set_xlabel(r'$M (h^{-1} M_\odot )$')
dndlnm_ax0.legend(loc = "best", numpoints = 1)
dndlnm_ax1.plot(MassesinhoverM , -sumandndlnM /psu.dndlnM(M = Masses, ps = psfrompk, cosmo = M000)*M000.h**3.0 )
#dndlnm_ax1.plot(MassesinhoverM , -sumanfsigma*sumandlsidlm*M000.Om0*sumanrhoc /Masses /psu.dndlnM(M = Masses, ps = psfrompk, M000 = M000) )
#dndlnm_ax1.plot(MassesinhoverM , -sumanfsigma*sumandlsidlm*M000.Om0*sumanrhoc /Masses /psu.dndlnM(M = Masses, ps = psfrompk, cosmo = M000) )
dndlnm_ax1.set_xscale('log')
plt.show()
sys.exit()
sys.exit()
plt.figure()
plt.plot(MassesinhoverM , massfnpkM000, label = "PS")
plt.plot(MassesinhoverM ,  sumandndlnM , label = "Suman")
plt.legend(loc="best")
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.figure()
plt.plot(MassesinhoverM , massfnpkM000/sumandndlnM )
plt.xlim(xmax = 1.e15)
plt.xscale('log')
plt.show()
sys.exit()
	#Obtain sigma from the power spectrum should be very close to 0.8 
#print "sigma ", psu.sigma(psfrompk, cosmo = Model)
	#Should be very close to 1.0 
#asrel  =  psu.getAsrel (psfrompk , sigma8 = 0.8, cosmo = Model)
#psfrompk2  = (psfrompk[0],asrel * psfrompk[1])
#print asrel
#Masses = np.logspace(10, 16, 40)
#plt.plot(Masses/0.71, psu.sigmaM(Masses/0.71, psfrompk,cosmo=Model)/psu.sigmaM(Masses/0.71, psrunfrompk,cosmo=Model), label =" run=-0.03")
#plt.legend(loc= "best",numpoints = 1)
#plt.xscale('log')
#plt.grid(True)
#plt.savefig("ratioofsigmaMfromPS.pdf")
#plt.show()
print "sigmaM ", psu.sigmaM(1.0e14/0.71, psfrompk2, cosmo = Model)
print "rhom", "{0:7.4e}".format(psu.__rhobg(cosmo = Model,unittype ="solarmassperMpc3")) , "units of solar mass/Mpc^3"
r = psu.filterradiusformass(1.0e14/0.71, cosmo = Model)
print "rad ", psu.filterradiusformass(1.0e14, cosmo = Model) , " in Mpc"
print "rad ", psu.filterradiusformass(1.0e14/0.71, cosmo = Model)*0.71 , " in Mpc/h for a mass of 1.0e14 in units of solar Mass/h"
print psu.__rhobg(cosmo= Model)
print psu.critdensity(unittype= "solarmassperMpc3", h = 0.71)
print "dlnsigmadlnR" , psu.dlnsigmadlnR (R = r,  ps = psfrompk2, z = 0.0 ,cosmo = Model) 
print "dRdlnM ", psu.dRdlnM (M = 1.0e14/0.71, cosmo = Model, z = 0.)
print "dlnsinv/dlnM", "{0:7.4e}".format(psu.dlnsigmainvdlnM ( M = 1.0e14 /0.71, ps = psfrompk, cosmo = Model))
print "dn/dlnM", "{0:7.4e}".format(psu.dndlnM ( M = 1.0e14 /0.71, ps = psfrompk, cosmo = Model))
sys.exit()
#print psu.dlnsigmainvdlnM ( M = 1.0e14 *0.71, ps = psfrompk, cosmo = Model)
import matplotlib.pyplot as plt
import numpy as np
R  = np.arange(1.0, 20.0,0.1)
#plt.plot(R, psu.sigma(ps = psfrompk2, R= R, cosmo = Model))
#plt.show()
print psu.__rhobg(cosmo=Model)/psu.critdensity(unittype= "solarmassperMpc3", h = 0.71)


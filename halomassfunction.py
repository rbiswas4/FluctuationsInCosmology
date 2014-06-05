#!/usr/bin/env python


#import utils
import numpy as np
import psutils  as psu
import massfunctions as mf


def ____dndlnM ( M ,
	ps , 
	z = 0. ,
	khmin = 1.0e-5,
	khmax = 2.0 ,
	logkhint = 0.005 ,
	bgtype = "cb", 
	powerspectrumfile = "LCDM_matterpower.dat" ,
	cosmo = None, 
	deltac = 1.686 , 
	fittingform = "Bhattacharya10",
	**params ):
	"""
	returns the mass function dn/dln(M) in units of h^3 Mpc^{-3}  
	args:
		M: mandatory, arraylike 
			mass bin in units of solar Mass
		powerspectrumfile : optional, string, defaults to 
			LCDM_matterpower.dat
			name of the power spectrum file from CAMB 
		cosmo: optional defaults to Planck13
			cosmology model

	returns:
		numpy array containing mass function in units of Mpc^{-3} 

	CHANGES:
		added argument deltac with default value 1.674 
		changed back to 1.686 R. Biswas 
	"""

	h = cosmo.H0/100.0
	
	sig = psu.sigmaM (M ,  
		ps , 
		bgtype = bgtype, 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint ,
		z = 0. ,
		cosmo = cosmo ,
		**params)

	sigm = sig

	dlsinvdlM =  -psu.dlnsigmadlnM (M ,
		ps ,
		z = 0., 
		bgtype = bgtype , 
		cosmo = cosmo ,
		khmin  = khmin ,
		khmax  = khmax ,
		logkhint = logkhint , 
		**params ) 

	if fittingform == "Bhattacharya10":
		f_sigma = mf.__fsigmaBhattacharya (
		 sigma = sigm,
		 deltac = deltac ,
		 z = z ,
		 A0 = 0.333 ,
		 a0 = 0.788 ,
		 p0 = 0.807 ,
		 q0 = 1.795 ,
		 alpha1 = 0.11 ,
		 alpha2 = 0.01 ,
		 alpha3 = 0.0 ,
		 alpha4 =  0.0,
		 Mlow = 6e11 ,
		 Mhigh = 3e15)    
	elif fittingform == "MICE" : 
		f_sigma = mf.fsigmaMICE(sigma = sigm, z  = z)
	else: 
		raise ValueError("This fitting form is not implemented")


	rhobg = psu.__rhobg( z =0. , bgtype = bgtype, 
		unittype = "solarmassperMpc3",  cosmo = cosmo)
	 
	#dndlnM = rhobg *sigm *dlsinvdlM /M 
	dndlnM = rhobg *f_sigma *dlsinvdlM /M 

	#dndlnM = dlsinvdlM  *f_sigma/M * rhobg
		#critdensity(h = cosmo.h, unittype = "solarmassperMpc3")*cosmo.Om0 
	#dndlnM =  rhobg * np.ones(len(sigm))
	#dndlnM =  rhobg * f_sigma
	return dndlnM   
def dndlnM ( M ,
	ps , 
	z = 0. ,
	khmin = 1.0e-5,
	khmax = 2.0 ,
	logkhint = 0.005 ,
	bgtype = "cb", 
	#powerspectrumfile = "LCDM_matterpower.dat" ,
	cosmo = None, 
	deltac = 1.674 , 
	fittingform = "Bhattacharya10",
	**params ):
	"""
	returns the mass function dn/dln(M) in units of h^3 Mpc^{-3}  
	args:
		M: mandatory, arraylike 
			mass in units of solar Mass
		ps: provide power spectrum tuple from CAMB at z = 0. The 
			required tuple is (koverh, ps)  
		cosmo: cosmology model of type FCPL 
		deltac : float 
			delta c value used in the fitting function
		fittingform: string, optional defaults to "Bhattacharya10"
			choices:
			"Bhattacharya10": Bhattacharya 2010 fitting form with 
				values suggested in paper"
			"MICE" 		: Crocce etal fitting form with 
				values suggested in paper
		

	returns:
		numpy array containing mass function in units of Mpc^{-3} 

	CHANGES:
		added argument deltac with default value 1.674 
	"""

	#z = np.asarray(z, dtype = float)

	h = cosmo.H0/100.0
	#rhocr = critdensity( h = h , 
	#	unittype = "solarmassperMpc3") 
	
	sig = psu.sigmaM (M ,  
		ps , 
		bgtype = bgtype, 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint ,
		z = 0.,
		cosmo = cosmo ,
		**params)

	sigm = sig
	#sigm = cosmo.growth(z=z)[0]*sig
	#print sigm

	dlsinvdlM =  -psu.dlnsigmadlnM (M ,
		ps ,
		z = 0. , 
		bgtype = bgtype , 
		cosmo = cosmo ,
		khmin  = khmin ,
		khmax  = khmax ,
		logkhint = logkhint , 
		**params ) 

	#used this to figure out what problems were
	#print z, sigm, dlsinvdlM , deltac 
	if fittingform == "Bhattacharya10":
		f_sigma = mf.__fsigmaBhattacharya (
		 sigma = sigm,
		 deltac = deltac ,
		 z = z ,
		 A0 = 0.333 ,
		 a0 = 0.788 ,
		 p0 = 0.807 ,
		 q0 = 1.795 ,
		 alpha1 = 0.11 ,
		 alpha2 = 0.01 ,
		 alpha3 = 0.0 ,
		 alpha4 =  0.0,
		 Mlow = 6e11 ,
		 Mhigh = 3e15)    
	elif fittingform == "MICE" : 
		f_sigma = mf.fsigmaMICE(sigma = sigm, z  = z)
	else: 
		raise ValueError("This fitting form is not implemented")

	rhobg = psu.__rhobg( z =0. , bgtype = bgtype, 
		unittype = "solarmassperMpc3",  cosmo = cosmo)
	 
	#dndlnM = rhobg *f_sigma *dlsinvdlM /M 
	#dndlnM = rhobg *sigm *dlsinvdlM /M 

	dndlnM = dlsinvdlM  *f_sigma/M * rhobg
		#critdensity(h = cosmo.h, unittype = "solarmassperMpc3")*cosmo.Om0 
	#return dndlnM 
	#dndlnM = rhobg*f_sigma 
	#dndlnM = rhobg*np.ones(len(sigm))
	#dndlnM = mf.__fsigmaBhattacharya(sigma = sigm, deltac = 1.686, z  = z)
	return dndlnM
def dndlnM0 ( M ,
	ps , 
	z = 0. ,
	khmin = 1.0e-5,
	khmax = 2.0 ,
	logkhint = 0.005 ,
	bgtype = "cb", 
	#powerspectrumfile = "LCDM_matterpower.dat" ,
	cosmo = None, 
	deltac = 1.674 , 
	fittingform = "Bhattacharya10",
	**params ):
	"""
	returns the mass function dn/dln(M) in units of h^3 Mpc^{-3}  
	args:
		M: mandatory, arraylike 
			mass in units of solar Mass
		ps: provide power spectrum tuple from CAMB at z = 0. The 
			required tuple is (koverh, ps)  
		cosmo: cosmology model of type FCPL 
		deltac : float 
			delta c value used in the fitting function
		fittingform: string, optional defaults to "Bhattacharya10"
			choices:
			"Bhattacharya10": Bhattacharya 2010 fitting form with 
				values suggested in paper"
			"MICE" 		: Crocce etal fitting form with 
				values suggested in paper
		

	returns:
		numpy array containing mass function in units of Mpc^{-3} 

	CHANGES:
		added argument deltac with default value 1.674 
	"""

	#z = np.asarray(z, dtype = float)

	h = cosmo.H0/100.0
	#rhocr = critdensity( h = h , 
	#	unittype = "solarmassperMpc3") 
	
	sig = psu.sigmaM (M ,  
		ps , 
		bgtype = bgtype, 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint ,
		z = 0.,
		cosmo = cosmo ,
		**params)

	#sigm = sig
	sigm = cosmo.growth(z=z)[0]*sig
	#print sigm

	dlsinvdlM =  -psu.dlnsigmadlnM (M ,
		ps ,
		z = 0. , 
		bgtype = bgtype , 
		cosmo = cosmo ,
		khmin  = khmin ,
		khmax  = khmax ,
		logkhint = logkhint , 
		**params ) 

	#used this to figure out what problems were
	#print z, sigm, dlsinvdlM , deltac 
	if fittingform == "Bhattacharya10":
		f_sigma = mf.__fsigmaBhattacharya (
		 sigma = sigm,
		 deltac = deltac ,
		 z = z ,
		 A0 = 0.333 ,
		 a0 = 0.788 ,
		 p0 = 0.807 ,
		 q0 = 1.795 ,
		 alpha1 = 0.11 ,
		 alpha2 = 0.01 ,
		 alpha3 = 0.0 ,
		 alpha4 =  0.0,
		 Mlow = 6e11 ,
		 Mhigh = 3e15)    
	elif fittingform == "MICE" : 
		f_sigma = mf.fsigmaMICE(sigma = sigm, z  = z)
	else: 
		raise ValueError("This fitting form is not implemented")

	rhobg = psu.__rhobg( z =0. , bgtype = bgtype, 
		unittype = "solarmassperMpc3",  cosmo = cosmo)
	 
	dndlnM = rhobg *f_sigma *dlsinvdlM /M 
	#dndlnM = rhobg *sigm *dlsinvdlM /M 

	#dndlnM = dlsinvdlM  *f_sigma/M * rhobg
		#critdensity(h = cosmo.h, unittype = "solarmassperMpc3")*cosmo.Om0 
	#dndlnM = M * dlsinvdlM
	#dndlnM = rhobg*np.ones(len(sigm))
	#dndlnM = mf.__fsigmaBhattacharya ( sigma= sigm, deltac = 1.686, z = z)

	return dndlnM   

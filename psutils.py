#!/usr/bin/env python
		
#CHANGES:
#Only assign values to cosmo as default if values from cosmo are being used
#otherwise pass as whole
#
#Fixed the spelling of filterradiusformass and its calls. Checked that the 
#tests at the bottom of the file still work. 
#R. Biswas, Thu Nov 14 15:31:46 CST 2013
#
#Fixed bug in function sigmaM, where filterradius is called without 
#optional argument cosmo and z. The same bug was there in sigmaM (where 
#it was being called without z, and derivatives of sigma and the mass 
# function calculation.
#R. Biswas, Thu Nov 14 17:59:14 CST 2013

import typeutils as tu
import massfunctions as mf
import numpy as np
import camb_utils.cambio as cio
import utils.filters as filters

def critdensity(h = 1.0, 
	unittype = 'kgperm3') :

	"""Returns the cosmology independent critical density as a function of 
	h = H0/100  ( 10^4 h^2 (3.0/(8.0 \pi G) )) in units of  kg /m^3. The
	answer is ~ 10^{-27} while the answer in gm/cm^3 is ~10^{-30} as quoted 
	in PDG for example. 	

	args:
		h: float, optional, defaults to 1.0 
		= H0/100 

		unittype: string ,optional, defaults to SI

			kgperm3: SI units 
			solarmassperMpc3: units of M_{\sun} / Mpc^3

	returns:
		critical density. If no argument is supplied the correct
		critical density is the return times h^2

	example usage:
		>>> from astropy.cosmology import Planck13 as cosmo
		>>> print critdensity (cosmo.h )
		>>> print cosmo.critical_density0
		>>> #The two above values should be same to unit definition 


	status:
		Tested with above tests and found to work. The value of 10^{-27}
		is due to the difference in units with astropy. 
		R. Biswas, Sun Aug 18 12:41:42 CDT 2013

		BUGFIX: unittype ought to be kgmerm3 in if loop, but was 
		written as kmperm3. Fixed Typo.
		R. Biswas,  Wed Nov 13 19:22:28 CST 2013

	"""

	from astropy import units as u
	from astropy import constants as ct 
	import numpy as np


	kmperMpc = (u.km / u.Mpc).decompose().scale 

	H0 = 100.0 * kmperMpc *h  # in units of per sec

	rhocrit = H0*H0 * 3.0 /(8.0 * np.pi * ct.G.value) 
	#Multiply mass in kg by convtosolmass to get mass in solar mass
	convtosolmass =  u.kg.to(u.solMass) 
	#Multiply distance in m to distance in Mpc
	convtoMpc = (u.m.to(u.Mpc))
	if unittype == "kgperm3" :
		rhocritu = rhocrit
	if unittype == "solarmassperMpc3":
		rhocritu = rhocrit*convtosolmass/convtoMpc**3.0 

	return rhocritu 

def __rhobg ( z = 0.0 , bgtype = "matter", unittype = "solarmassperMpc3", 
	cosmo = None):

	"""returns the background density at redshift z. If bgtype = "matter"
	then the function returns the background matter (CDM + Baron + massive 
	neutrino)  density. 

	args:
		z: 
			optional, float , defaults to 0.
			redshift at which the background density is 
			calculated

		bgtype :
			string, optional, defaults to "matter"
			choices: "matter"
				"cb" : baryon - cdm 
		
		unittype : 
			string, optional ,defaults to "Msun/Mpc3"
			defines the unit of the return 

			solarmassperMpc3:
				Units of solar mass per Mpc cube

			SI :
		cosmo	: w0wa cosmology
			
	returns:
		float, the background density in units of type unittype

	example usage:
	status:
	notes: 
		Will change later. This exists because our definition of 
		matter includes neutrinos, which I want to remove later 
		so that matter is baryon CDM by definition. 
		
	""" 

	#if cosmo == None:
	#	from astropy.cosmology import Planck13 as cosmo

	if bgtype == "matter":
		Om0 = cosmo.Om0
	if bgtype == "cb":
		Om0 = cosmo.Ob0 + cosmo.Oc0

	h = cosmo.H0 /100.
	rho =  critdensity(h = h, unittype = 'solarmassperMpc3')*Om0*(1.+z)**3.0
	return rho

def filterradiusformass( M , 
	z = 0. , 
	bgtype = "matter" , 
	cosmo = None):

	"""
	returns a radius in units of Mpc which encloses a mass M of the
	homegenous density of particles specified as bgtype at the redshift 
	z for a cosmology cosmo.   
	args:
		M : 
			mass in solar masses  
			z : float, defaults to 0. 
			redshift 
		bgtype  : string, defaults to matter 
			background density to use in converting 
			mass to a radius
	cosmo : wowa cosmology, defaults to Planck13
	returns :
		Radius in Mpc
	"""
	import numpy as np

	#if cosmo == None:
	#	from astropy.cosmology import Planck13 as cosmo
		
	rhobg =  __rhobg ( z , bgtype = bgtype, unittype = "solarmassperMpc3", cosmo = cosmo)

		#Assume mass is made of all matter
	Rcube = 3* M / rhobg /4.0 / np.pi 
	R = np.power (Rcube, 1.0/3.0)

	return R
	
def powerspectrumfromfile(fname, 
	koverh = None ,
	pstype = "matter" , 
	h  = 0.71 ,
	ns = 0.963 ,
	As = 2.1e-9 , 
	Omegacdm = None , 
	Omegab   = None ):

	"""
	returns a tuple of koverh values and the interpolated power
	spectra at these values of hoverh using a CAMB output which 
	may be a power spectrum output or transfer function output. 
	If the output is a transfer function output, then ns,  h, 
	and As must be supplied 

	args:

	returns:
		tuple of (koverh, ps values)

	"""

	#decide if file is transfer function or Power spectrum output

	psfile = False
	tffile = False
	Unknown = True

	tmpfile = np.loadtxt(fname)
	shapetuple = np.shape(tmpfile)
	if shapetuple[-1] == 7:
		tffile  = True
		Unknown = False 
	if shapetuple[-1] ==2 :
		psfile  = True 	
		Unknown  = False
	if koverh == None:
		koverh = tmpfile[:,0]
	

	if Unknown:
		#file is not CAMB transfer function or power spectrum output
		raise ValueError("Unknown filename supplied")

	if psfile :
		if pstype != "matter" :
			raise ValueError ("Cannot obtain non-matter power spectrum  from CAMB power spectrum file")
		return (koverh , powerspectrum(koverh, fname )	)


	if tffile :
		if pstype == "matter" :
			transfer  = cio.loadtransfers(rootname = None, 
				filename = fname) 
			ps = cio.matterpowerfromtransfersforsinglespecies( 				koverh ,
				transfer , 
				h ,
				As ,
				ns )


			return (ps [:,0], ps[:,1])

		elif pstype == "cb"    :

			#LOST CODE HERE
			return 0
def powerspectrum ( koverh , powerspectrumfile = "LCDM_matterpower.dat" ):

	"""
	returns interpolated values of the powerspectrum in the 
	powerspectrumfile with k values in units of h/Mpc
		P(k) is in 
	NOTS URE
	"""

	import numpy as np

	pk = cio.loadpowerspectrum(powerspectrumfile)

	return  np.interp( koverh, pk[:,0],pk[:,1])



def sigma(ps ,  R = 8 ,  khmin = 1e-5, khmax = 2.0, logkhint = 0.005, cosmo = None, filt = filters.Wtophatkspacesq, **params) :



	"""
	usage:
		>>> pk = np.loadtxt("powerspectrum") 
		>>> sigma (ps = (pk[:,0],pk[:,1]), cosmo = cosmo)

	"""

	sigsq= sigmasq(ps = ps , R  =R, khmin = khmin , khmax = khmax , 
		logkhint = logkhint , cosmo=cosmo , filt = filt , **params )

	return np.sqrt(sigsq )
def sigmasq (ps , R = 8. , khmin = 1.0e-5 , khmax = 2.0, logkhint = 0.005 , 
	cosmo = None, filt= filters.Wtophatkspacesq,  **params) : 
	"""
	Returns the variance of the overdensity field smoothed at 
	a radius of R Mpc/h using a filter specified by filt

	args:
		ps: tuple of koverh, power spectrum values


	notes:
		We need h, even if CAMB power spectrum is given
	"""

	import numpy as np
	import scipy.integrate as si

	logkhmin  = np.log(khmin)
	logkhmax  = np.log(khmax)
	logkh = np.arange(logkhmin, logkhmax , logkhint) 

	h = cosmo.H0/100.0
	khvals = np.exp(logkh)
	k  = khvals * h  
	

	psinterp = np.interp (khvals , ps[0], ps[1])

	#print khvals
	#print psinterp
	if tu.isiterable(R):
		R = np.asarray(R)
		kr = np.outer( R, khvals ) 
	else:
		kr  = R* khvals 
	kwinsq=   filt (kr, R)
	#kwin =   W (kr)
	ksqWsqPk = k*k *kwinsq* psinterp /2. /np.pi/ np.pi/h /h/h	
	
	
	sigmasq = si.simps ( ksqWsqPk, x = k) 
	return sigmasq	


def getAsrel (ps , sigma8, khmin = 1.0e-5 , khmax = 2.0, logkhint = 0.005 , 
	cosmo = None, filt= filters.Wtophatkspacesq,  **params) :
	"""
	returns a relative value of As by which to multiply the power spectrum
	values in order to obtain sigma8
	args:

	returns:
		float, Asrel
	"""

	sigsq = sigmasq (ps , khmin= khmin, khmax =khmax,logkhint =logkhint, cosmo = cosmo, filt = filt , **params)  	
	Assq = sigma8*sigma8 / sigsq 
	return np.sqrt(Assq)
 

def sigmaM (M , 
	ps , 
	bgtype = "matter", 
	khmin = 1.0e-5 , 
	khmax = 2.0  ,
	logkhint = 0.005 ,
	z = 0.0 ,
	cosmo = None ,
	**params):

	"""Returns the standard deviation of the overdensity fields 
	smoothed at a radius corresponding to mass M.
	args:
		M: array like , mandatory
			mass of halo in units of solar masses
		
	notes:
		the bgtype matters for converting Mass M to R.
		Also ps must be set accordingly

	"""	

	if tu.isiterable(M):
		M =  np.asarray(M)

	#if cosmo == None:
	#	from astropy.cosmology import Planck13 as cosmo

	h = cosmo.H0/100.0

	R =  filterradiusformass( M  , bgtype= bgtype, z = z, cosmo = cosmo)
	RinMpcoverh = R*h 

	return sigma( ps , R = RinMpcoverh, khmin = khmin , khmax = khmax, logkhint = logkhint , cosmo= cosmo, **params) 


		
def dlnsigmainvdlnM (M ,
	ps ,
	z = 0.0 , 
	bgtype = "matter", 
	cosmo = None , 
	khmin  = 1.0e-5 ,
	khmax  = 2.0 ,
	loghkint = 0.005 , 
	**params ) :	
	"""
	returns the derivative dln (\sigma^{-1})/ d ln M at values of M by

	args:
		M: array-like, mandatory 
			mass of halo in units of solar mass
		ps : tuple, mandatory
			(koverh , ps)

	notes:
		dln(\sigma^{-1})/dln M = M /sigma^{-1}* (-1.)\sigma^{-2}* d sigma /dR dR/dM 
		Done by calculating 1/sigma * dsigma /dR * dR /d ln M  , 
		where dsigma / dR = sigma with derivative of filter  

	"""

	sig = sigmaM (M ,  ps , khmin = 1.0e-5 , 
		khmax = 2.0  ,
		logkhint = 0.005 ,
		z = 0.0 ,
		bgtype = bgtype , 
		cosmo = cosmo ,
		**params)

	h = cosmo.h
	R =  filterradiusformass( M  , bgtype = bgtype, z = z, cosmo = cosmo)
	dRdlnM = M * (1.0/3.0)/R 
	RinMpcoverh = R*h 

	dsigdR = dsigmadR( R= RinMpcoverh, z = z , ps = ps, bgtype= bgtype, cosmo = cosmo, khmin = khmin, khmax = khmax, loghkint = loghkint)

	return -1.0/ sig * dsigdR * dRdlnM  

def dRdM ( M , 
	bgtype = "matter", 
	cosmo = None ):
	"""
	args: 
		M : 
		in units of solar masses
		

	"""

	R = filterradiusformass( M  , bgtype = bgtype, z = z, cosmo = cosmo)

	dRdM = (1.0/3.0)/R 

	return dRdM 


def dsigmadR (R ,
	ps ,
	z = 0.0 , 
	bgtype = "matter", 
	cosmo = None , 
	khmin  = 1.0e-5 ,
	khmax  = 2.0 ,
	loghkint = 0.005 , 
	**params ) :	
	"""
	returns the derivative d(\sigma)/ dR  at values of M by

	args:
		M: array-like, mandatory 
			mass of halo in units of solar mass
		ps : tuple, mandatory
			(koverh , ps)

	notes:
		Done by calculating 1/sigma * dsigma /dR * dR /d ln M  , 
		where dsigma / dR = sigma with derivative of filter  

	"""

	dsigmasqdR = sigmasq(ps ,  R = R ,  khmin = 1e-5, khmax = 2.0, logkhint = 0.005, cosmo = cosmo, filt = filters.dWtophatkspacesqdR , **params) 
	sig = sigma(ps ,  R = R ,  khmin = 1e-5, khmax = 2.0, logkhint = 0.005, cosmo = cosmo, filt = filters.Wtophatkspacesq , **params) 
	

	return dsigmasqdR/2.0 / sig 

def dlninvsigmaMdlnM ( M ,
	ps , 
	khmin = 1.0e-5,
	khmax = 2.0 ,
	z = 0.0 , 
	logkhint = 0.005 ,
	cosmo = None, 
	**params ):

	"""
	returns the derivative dln(sigma^{-1}) / d lnM at values 
	of M except the edge values.   
	args:
		M: numpy array of mass:
	notes:
		Expect M to be equally spaced in log space
	"""

		#expect points to be equally 
		#spaced in logspace
	lnM = np.log(np.asarray( M))
		#Displace points to middle of the the input 
		#array for evaluations	
	evalpoints = 0.5*(lnM[:-1] + lnM[1:])
		#Exponentiate back to mass values , rather than log mass
	Meval = np.exp(evalpoints)

	#import matplotlib.pyplot as plt

	#plt.plot (Meval , M[1:])
	#plt.show()
	
		#Compute inverse sigma values for the masses
	invsigmavals = 1./sigmaM( Meval , 
		z = z , 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint  ,
		powerspectrumfile = powerspectrumfile ,
		cosmo = cosmo ,
		**params)

		# Get ln 
	loginvsigmavals = np.log(invsigmavals) 
		#Differences and derivatives
	deltaloginvsigmavals = loginvsigmavals[1:] - loginvsigmavals [:-1]
	deltalogM = evalpoints[1:] - evalpoints[:-1] 

	return deltaloginvsigmavals / deltalogM 
	

def dndlnM ( M ,
	z = 0. ,
	khmin = 1.0e-5,
	khmax = 2.0 ,
	logkhint = 0.005 ,
	powerspectrumfile = "LCDM_matterpower.dat" ,
	cosmo = None, 
	deltac = 1.674 , 
	**params ):
	"""
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
	"""

	if cosmo == None:
		from astropy.cosmology import Planck13 as cosmo

	h = cosmo.H0/100.0
	rhocr = critdensity( h = h , 
		unittype = "solarmassperMpc3") 

	dlninvsigmadlnMvals = dlninvsigmaMdlnM ( M ,
		z = z , 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint  ,
		powerspectrumfile = powerspectrumfile ,
		cosmo = cosmo ,
		**params)
		
	sigmaMval  = sigmaM (M[1:-1], 
		z = z , 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint  ,
		powerspectrumfile = powerspectrumfile ,
		cosmo = cosmo ,
		**params)

	f_sigma = mf.__fsigmaBhattacharya (
		 sigma = sigmaMval,
		 deltac = deltac ,
		 z = 0.0 ,
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



	rhobg = __rhobg( z =z , bgtype = "matter", 
		unittype = "solarmassperMpc3",  cosmo = cosmo)
	dndlnM = rhobg *f_sigma *dlninvsigmadlnMvals /M[1:-1]

	return dndlnM   
if __name__=="__main__":

	import numpy as np
	import matplotlib.pyplot as plt
	import camb_utils.cambio  as cio
	import sys
	#pk = cio.loadpowerspectrum ("example_data/LCDM_def_matterpower.dat")
	pk = cio.loadpowerspectrum ("LCDM_matterpower.dat")
	ts = cio.loadtransfers(filename = "example_data/LCDM_def_transfer_out.dat")
	print np.shape(ts)
	#print pk[:,0]
	pkt = cio.matterpowerfromtransfersforsinglespecies(koverh = pk[:,0], 
		transfer = (ts[:,0],ts[:,-1]), h = 0.71, As = 2.1e-9, ns  = 0.963)
	
	plt.loglog ( pk[:,0], pk[:,1])	
	plt.loglog ( pkt[:,0], pkt[:,1])
	
	plt.figure()
	from astropy.cosmology import Planck13 as cosmo
	print sigma(ps = (pk[:,0],pk[:,1]) , R = 8.0, cosmo = cosmo)
	#plt.show()
	sys.exit()
	M = 10.**(np.arange(7,16,0.2))
	R = np.arange(0.0005, 50.0,0.1)
	#R = np.array([4,8,12])
	print sigma (8.0)
	plt.figure()
	plt.plot(R, sigma(R))
	plt.xlabel("R ( Mpc /h )") 
	plt.ylabel(r'$\sigma (R)$') 
	plt.figure()
	plt.plot(M ,filterradiusformass(M))
	plt.xscale('log')
	plt.xlabel("M ") 
	plt.ylabel(r'$R(M) Mpc $') 
	plt.figure()
	plt.plot ( M, sigmaM(M, powerspectrumfile = "LCDM_def_matterpower.dat"), "o")
	plt.plot ( M, 1./sigmaM(M,powerspectrumfile = "LCDM_def_matterpower.dat"),'o')
	plt.xlabel("M ") 
	plt.ylabel(r'$\sigma (M)$') 
	plt.xscale('log')
	plt.figure()
	plt.plot ( M[1:-1], dlninvsigmaMdlnM (M ),"o") 
	plt.xlabel(r'$M (M_\odot$') 
	plt.ylabel(r'$\frac{d ln \sigma^{-1}}{d ln(M)}$') 
	plt.xscale('log')
	plt.tight_layout()

	plt.savefig("dlninvsigmadlnM,pdf")

	plt.figure()
	#plt.plot (1./ sigmaM(M[1:-1]), dndlnM (M ), "o")
	plt.plot (M[1:-1], dndlnM (M ), "o")
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	#print filterradiusformass ( M = 
	#plt.show()

	

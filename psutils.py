#!/usr/bin/env python
#
#This is a set of wrappers designed to use methods of obtaining linear 
#quantities of interest from outputs of actual programs taht do the 
#calculations, like CAMB with the help of utilities for specific programs.
#
#USEFUL ROUTINES:
#
#powerspectrum: obtains the linear power spectrum of various quantities from
#---------------
#	standard outputs of programs like CAMB
#sigma
#---------------
#sigmaM
#--------------
		
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

import sys
import matplotlib.pyplot as plt
import utils.typeutils as tu
import massfunctions as mf
import growthfunction
import numpy as np
import camb_utils.cambio as cio
import utils.filters as filters


verbose  = False
def loginterp(xreq, xnative, ynative , left = np.nan , right = np.nan):

        logxreq = np.log(xreq)
        npinterp  = np.interp(logxreq , np.log(xnative), np.log(ynative), left = np.nan, right = np.nan)

        return np.exp(npinterp)

def critdensity(h = 1.0, 
	unittype = 'kgperm3') :

	"""Returns the critical density today $\Omega_{crit}(z=0) as a 
	function of h through the formula 
	$10^4 h^2 (3.0/(8.0 \pi G) )) in units of  kg /m^3 or solar 
	masses / Mpc^3 . 

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
	notes : 
		The answer is ~ 10^{-27} while the answer in gm/cm^3 is 
		~10^{-30} as quoted in PDG for example. 	
		
		TODO: This function will be removed to cosmodefs

	"""

	from astropy import units as u
	from astropy import constants as ct 


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

def __rhobg ( z = 0.0 , bgtype = "cb", unittype = "solarmassperMpc3", 
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

		Remove to cosmodefs. 
		
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
	bgtype = "cb" , 
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
	****************************************
	DEPRECATED: USE cambio functions instead
	*****************************************
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
def powerspectrum ( koverh , 
	asciifile = None ,
	pstype = "matter", 
	sigma8type = "matter" ,
	method = "CAMBoutfile",
	z      = 0.0 , 
	cosmo =  None ,
	interpmethod  = 'log' , 
	**params):

	"""
	returns linearly interpolated values of the powerspectrum in the 
	powerspectrumfile with k values in units of h/Mpc. Using
	this with koverh = None, returns the values in the table. 

	args:
		koverh : array-like of floats or Nonetype, mandatory
			k in units of h/Mpc
		asciifile: string, 
			Filename for power spectrum or CAMB transfer function. 
			power sepctrum or transfer function input will be 
			recognized from CAMB file structure.
		cosmo    : interfacecosmology/astropy cosmological model
			
		method   : string, optional , defaults to "CAMBoutfile" 
			Method of obtaining power spectrum with fixed options
			options:
			-------
			CAMBoutfile   :assume that the asciifile output of CAMB 
				is at desired redshift 
			CAMBoutgrowth :Use the asciifile from CAMB output at 
				z  = 0.0 , and use a growth function to find 
				the power spectrum at z = z
		interpmethod: string, optional, defaults to 'log'
			options:
				'log': linearly interpolates 
				log(koverh) vs log(PS) in units of h^#/Mpc^3
				'linear' : linearly interpolates koverh and PS
	returns:
		tuple (koverh , power spectrum)

	notes: should be able to obtain the powerspectrum in a variety of 
		methods with code being added
		Override Rules: 
			sigma8 overrides As 
			params dictionary overrides cosmo 
	"""

	#Make sure evaluation method is implemented
	if not method in ["CAMBoutfile","CAMBoutgrowth"]:
		raise ValueError("Method not defined")
	if method in ["CAMBoutfile", "CAMBoutgrowth"] :
		#Query CAMB file type: power spectrum or transfer
		psfile, tkfile, Unknown = cambasciifiletype (asciifile)  

	if method == "CAMBoutgrowth":
		#cannot calculate growth if cosmo not provided
		if cosmo == None and z >0.000001 :
			raise ValueErrror("Method does not work if cosmo is not defined")

	Dg = 1.0 
	if cosmo !=None:
		if z > 0.000001:
			Dg = cosmo.growth(z)[0]
		
	# decide whether As or sigma8 is to be used
	# if sigma8 provided use it, otherwise As
	sigma8 = None
	As     = None 

	if params.has_key("sigma8"):
		if params.sigma8 != None : 
			sigma8 = params["sigma8"]
		if params.As != None :
			As = params["As"]

	if cosmo != None :
		if cosmo.sigma8  != None:
			if sigma8 == None: 
				sigma8 = cosmo.sigma8 
		if cosmo.As      != None:
			if As == None :
				As = cosmo.As

	#If neither As or sigma8 are provided fail!
	if As == None and sigma8 == None and not psfile :
		raise ValueError("without As or sigma8 provided, matter power spectrum cannot be calculated from transfer functions\n") 
	if sigma8 != None:
		As = 1.0
	

	if params != None:
		paramdict = params
	paramdict["As"] = As
		
	#print "VALUES passed on from powerspectrum routine \n"		
	if verbose:
		print "sigma8 = ", sigma8, " As ", As
	#print paramdict["As"], "IN powerspectrum"

	pstmp  = __powerspectrum ( koverh = None, 
		asciifile = asciifile  ,
		pstype = pstype ,
		method = method ,
		z      = z , 
		cosmo  =  cosmo ,
		**paramdict )

	#If sigma8 is given, we need to normalize power spectrum
	#power spectrum to normalize is pssigma8
	if sigma8 != None: 
		if pstype !=sigma8type:	
			if verbose:
				print "evaluate sigmatype ps \n"
			pssigma8  = __powerspectrum ( koverh = None, 
				asciifile = asciifile  ,
				pstype = sigma8type ,
				method = method ,
				z      = z , 
				cosmo  =  cosmo ,
				**paramdict)
	
		else:
			pssigma8  = pstmp 

	
	if sigma8 != None :

		Asrel =  getAsrel (pssigma8 , sigma8, cosmo = cosmo, 
			filt= filters.Wtophatkspacesq,  **paramdict) 
		#print "Now As has been determined to be ", sigma8type , Asrel
		v = pstmp[0], Asrel*pstmp[1] 
	else :
		v = pstmp 

	
	if koverh != None:
		if interpmethod == "linear":
			ret = koverh, np.interp(koverh, v[0], v[1], 
				left = np.nan , right = np.nan) 
		else:
			interpmethod = "log"
			ret = koverh, loginterp(koverh, v[0], v[1], 
				left = np.nan , right = np.nan) 
			
	else: 
		ret =  v 
 
	if method == "CAMBoutgrowth" :
		return ret[0],Dg*Dg*ret[1]
	else: 
		return ret
def getvalsfromparams(cosmo, **params):

 
	""" 
	TO DO 
	provide a general function to pass values into cosmo and params
	"""

	return None

def cambasciifiletype( fname ) :

	# Decide whether this ia a matter or transfer file
	psfile = False
	tkfile = False
	Unknown = True

	tmpfile = np.loadtxt(fname )
	shapetuple = np.shape(tmpfile)
	if shapetuple[-1] == 7:
		tkfile  = True
		Unknown = False 
	if shapetuple[-1] ==2 :
		psfile  = True 	
		Unknown  = False

	if Unknown:
		#file is not CAMB transfer function or power spectrum output
		raise ValueError("Unknown filename supplied")
		

	return psfile, tkfile, Unknown 

	
def __powerspectrum ( koverh , 
	asciifile = None ,
	pstype = "matter", 
	method = "CAMBoutfile",
	z      = 0.0 , 
	cosmo =  None ,
	**params):

	"""
	DO NOT CALL DIRECTLY. CALL powerspectrum instead 
	returns linearly interpolated values of the powerspectrum in the 
	powerspectrumfile with k values in units of h/Mpc. Using
	this with koverh = None, returns the values in the table. 

	args:
		koverh : array-like of floats or Nonetype, mandatory
			k in units of h/Mpc
		asciifile: string, 
			Filename for power spectrum or CAMB transfer function. 
			power sepctrum or transfer function input will be 
			recognized from CAMB file structure.
		method   : string, optional , defaults to "CAMBoutfile" 
			Method of obtaining power spectrum with fixed options
			options:
			-------
			CAMBoutfile   :assume that the asciifile output of CAMB 
				is at desired redshift 
			CAMBoutgrowth :Use the asciifile from CAMB output at 
				z  = 0.0 , and use a growth function to find 
				the power spectrum at z = z
			
			
	returns:
		tuple (koverh , power spectrum)

	notes: should be able to obtain the powerspectrum in a variety of 
		methods with code being added
	"""

	#ensure we are supposed to read CAMB outfiles
	if not method in ["CAMBoutfile","CAMBoutgrowth"]:
		raise ValueError("Method not defined")

#	# Decide whether this ia a matter or transfer file
	#This has been made a function
#	psfile = False
#	tkfile = False
#	Unknown = True
#
#	shapetuple = np.shape(tmpfile)
#	if shapetuple[-1] == 7:
#		tkfile  = True
#		Unknown = False 
#	if shapetuple[-1] ==2 :
#		psfile  = True 	
#		Unknown  = False


	psfile, tkfile, Unknown = cambasciifiletype ( asciifile )
	
	tmpfile = np.loadtxt(asciifile)
	if koverh == None:
		koverh = tmpfile[:,0]

	if Unknown:
		#file is not CAMB transfer function or power spectrum output
		raise ValueError("Unknown filename supplied")

	if psfile: 
		pk = cio.loadpowerspectrum(asciifile)
		if not np.all(np.diff(pk[:,0])>0.):
			raise ValueError("The k values in the power spectrum file are not in ascending order")

		if koverh == None :
			return (pk[:,0], pk[:,1])

		return  koverh, np.interp( koverh, pk[:,0],pk[:,1],left = np.nan, right = np.nan)
	if tkfile:
		#print "AS " , params["As"]
		#print cosmo.Ob0, cosmo.Oc0
	
		if pstype == "cb":
			#print "filename ", asciifile
			pk = cio.cbpowerspectrum ( transferfile = asciifile , 
				Omegacdm  =  cosmo.Oc0,
				Omegab    =  cosmo.Ob0, 
				h         =  cosmo.h, 
				Omeganu   =  cosmo.On0,
				As        =  params["As"], 
				#As        =  cosmo.As, 
				ns        =  cosmo.ns, 
				koverh    =  None )

			return (pk [:,0], pk[:,1])

		if pstype == "matter" :
			if koverh == None :
				koverh = tmpfile[:,0]

			transfer  = cio.loadtransfers( filename = asciifile) 

			transfertuple = (transfer[:,0], transfer[:,-1])
			ps = cio.matterpowerfromtransfersforsinglespecies(
				koverh ,
				transfer  = transfertuple,
				h = cosmo.h ,
				As = params["As"],
				ns = cosmo.ns)


			return (ps [:,0], ps[:,1])

		
		return koverh, pk 






def sigma(ps ,  R = 8 ,  khmin = 1e-5, khmax = 2.0, logkhint = 0.005, cosmo = None, filt = filters.Wtophatkspacesq, **params) :



	"""
	returns the square root of the variance of isotropic, homogeneous 
	fluctuations filtered with a single scale filter at a scale of 
	R Mpc/h. 
	args: 
		ps: tuple of koverh , power spectrum values 
		R : array-like float, optional defaults to 8.0
			radius in units of Mpc/h over which the filtering
			is done 
		filt: function describing the shape of the filter
			default is filters.Wtophatkspacesq which is 
			the Fourier transform of the tophat at R Mpc/h
		cosmo: cosmological model
	usage:
		>>> pk = np.loadtxt("powerspectrum") 
		>>> sigma (ps = (pk[:,0],pk[:,1]), cosmo = cosmo)

	"""

	sigsq= sigmasq(ps = ps , R  =R, khmin = khmin , khmax = khmax , 
		logkhint = logkhint , cosmo=cosmo , filt = filt , **params )

	return np.sqrt(sigsq )
def sigmasq (ps , R = 8. , usenative = True, khmin = 0.9e-5 , khmax = 5.0, logkhint = 0.005 , 
	cosmo = None, filt= filters.Wtophatkspacesq,  **params) : 
	"""
	Returns the variance of the overdensity field smoothed at 
	a radius of R Mpc/h using a filter specified by filt

	args:
		ps: tuple of koverh, power spectrum values
		R : float array like
			distance scale in units of Mpc/h over which 
			the filtering is done
		usenative: bool, optional , defaults to True
			Use values provided in ps, rather than 
			interpolation
		cosmo: Model, whose hubble constant will be used
		
	returns :
		array of sigmasq values 

	notes:
		- We need h, even if CAMB power spectrum is given
		- If interpolation is used only, and the range provided
			is outside the range of the data, only those points
			in the original range will be used. extrapolation
			is dangerous, particularly at high k, unless it is 
			made to drop as a power law. 

	"""

	import numpy as np
	import scipy.integrate as si


	h = cosmo.H0/100.0
	if usenative :
		khvals = ps[0]
	else: 
		logkhmin  = np.log(khmin)
		logkhmax  = np.log(khmax)
		logkh = np.arange(logkhmin, logkhmax , logkhint) 
		khvals = np.exp(logkh)

		logkhmin = max(min(ps[0]),logkhmin)
		logkhmax = min(max(ps[0]),logkhmax) 

	k  = khvals * h  
	

	psinterp = np.interp (khvals , ps[0], ps[1], left = np.nan, right = np.nan)

	#plt.loglog(khvals, psinterp, label="interp")
	#plt.loglog(ps[0], ps[1], label="native")
	#plt.legend(loc= "best")
	#plt.show()

	if tu.isiterable(R):
		R = np.asarray(R)
		kr = np.outer( R, khvals ) 
	else:
		kr  = R* khvals 
	kwinsq=   filt (kr, R)
	#kwin = 3*(np.sin(kr)-kr*np.cos(kr))/(kr)**3
	#kwinsq = kwin *kwin 
	
	ksqWsqPk = k*k *kwinsq* psinterp /2. /np.pi/ np.pi/h /h/h	
	
	
	sigmasq = si.simps ( ksqWsqPk, x = k, even = 'avg') 
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
	Asrel = sigma8*sigma8 / sigsq 
	return Asrel
 

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

	
	print "RinMpcoverh ***************"
	print RinMpcoverh 
	#return RinMpcoverh 	
	return sigma( ps , R = RinMpcoverh, khmin = khmin , khmax = khmax, logkhint = logkhint , cosmo= cosmo, **params) 


		
def dlnsigmadlnM (M ,
	ps ,
	z = 0.0 , 
	bgtype = "matter", 
	cosmo = None , 
	khmin  = 1.0e-5 ,
	khmax  = 2.0 ,
	logkhint = 0.005 , 
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

	sig = sigmaM (M ,  ps , khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint ,
		z = z ,
		bgtype = bgtype , 
		cosmo = cosmo ,
		**params)

	h = cosmo.h
	R =  filterradiusformass( M  , bgtype = bgtype, z = z, cosmo = cosmo)

	dlnRdlnM = 1.0/3.0
	RinMpcoverh = R*h 
	
	#d ln sigma /d ln R = d ln sigma^2 / d ln R / sigma^2/ 2.0  
		#sigmasq with filter of dWtophatkspacesqdlnR  
		# is dln sigma^2/ d ln R 
	dlnsigdlnR =  sigmasq (R = RinMpcoverh , ps  = ps, z = z ,  
		bgtype = bgtype,  filt = filters.dWtophatkspacesqdlnR, cosmo = cosmo ,  khmin  = khmin , 
		khmax  = khmax , logkhint = logkhint,  **params )/sig/ sig/2.0 

	#return sig
	return  dlnsigdlnR *dlnRdlnM


def dndlnM ( M ,
	ps , 
	z = 0. ,
	khmin = 1.0e-5,
	khmax = 2.0 ,
	logkhint = 0.005 ,
	bgtype = "matter", 
	powerspectrumfile = "LCDM_matterpower.dat" ,
	cosmo = None, 
	deltac = 1.674 , 
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
	"""

	h = cosmo.H0/100.0
	#rhocr = critdensity( h = h , 
	#	unittype = "solarmassperMpc3") 

	
	sig = sigmaM (M ,  
		ps , 
		bgtype = bgtype, 
		khmin = khmin , 
		khmax = khmax  ,
		logkhint = logkhint ,
		z = z,
		cosmo = cosmo ,
		**params)

	dlsinvdlM =  -dlnsigmadlnM (M ,
		ps ,
		z = z , 
		bgtype = bgtype , 
		cosmo = cosmo ,
		khmin  = khmin ,
		khmax  = khmax ,
		logkhint = logkhint , 
		**params ) 

	f_sigma = mf.__fsigmaBhattacharya (
		 sigma = sig,
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

	rhobg = __rhobg( z =z , bgtype = bgtype, 
		unittype = "solarmassperMpc3",  cosmo = cosmo)
	 
	dndlnM = rhobg *f_sigma *dlsinvdlM /M 

	#dndlnM = dlsinvdlM  *f_sigma/M * rhobg
		#critdensity(h = cosmo.h, unittype = "solarmassperMpc3")*cosmo.Om0 
	return dndlnM   
if __name__=="__main__":

	import numpy as np
	import matplotlib.pyplot as plt
	import camb_utils.cambio  as cio
	import sys
	#pk = cio.loadpowerspectrum ("example_data/LCDM_def_matterpower.dat")
	pk = cio.loadpowerspectrum ("LCDM_matterpower.dat")
	ts = cio.loadtransfers(filename = "example_data/LCDM_def_transfer_out.dat")
	#print np.shape(ts)
	#print pk[:,0]
	pkt = cio.matterpowerfromtransfersforsinglespecies(koverh = pk[:,0], 
		transfer = (ts[:,0],ts[:,-1]), h = 0.71, As = 2.1e-9, ns  = 0.963)
	
	plt.loglog ( pk[:,0], pk[:,1])	
	plt.loglog ( pkt[:,0], pkt[:,1])
	
	plt.figure()
	from astropy.cosmology import Planck13 as cosmo
	#print sigma(ps = (pk[:,0],pk[:,1]) , R = 8.0, cosmo = cosmo)
	#plt.show()
	sys.exit()
	M = 10.**(np.arange(7,16,0.2))
	R = np.arange(0.0005, 50.0,0.1)
	#R = np.array([4,8,12])
	#print sigma (8.0)
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

	

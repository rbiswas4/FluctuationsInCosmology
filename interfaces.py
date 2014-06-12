#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology  import w0waCDM 

def _hacccosmologydict(haccinputdict ) :

	cosmologyvariables =  ['HUBBLE', 'OMEGA_CDM','DEUT','OMEGA_NU' ,'W_DE','WA_DE', 'SS8', 'NS']
	hacccosmologydict = {}
	for key in haccinputdict :
		if key in cosmologyvariables :
			hacccosmologydict[key] = float(haccinputdict[key])

	return hacccosmologydict 
class FCPL (w0waCDM) :

	def __init__ ( self, 
		H0 = 70. , 
		Om0 = 0.3, 
		Ob0 = 0.02256, 
		Ok0 = 0.0  ,
		ns = 0.96 , 
		As = None , 
		sigma8 = None , 
		w0 = -1., 
		wa = 0., 
		sigmamnu = 0.0  ,
		Tcmb0=2.725,
                Neff=3.04, 
		TFfilelist  = [] , 
		TFredshiftlist = [],
		setvaluesfrom = None , 
		name='FCPL'):
		

		from astropy.cosmology import FlatwCDM 
		tmpmodel = FlatwCDM(
			H0 = H0, 
			Om0 = Om0, 
			#Ok0 = 0.,
			w0 = w0 , 
			Tcmb0 = Tcmb0,
			Neff =Neff, 
			name = 'tmp') 

		Ode0= 1.0 - Om0 - tmpmodel.Ogamma0 - tmpmodel.Onu0 - Ok0 

		w0waCDM.__init__(self, 
			H0 = H0, 
			Om0 = Om0, 
			Ode0 = Ode0 ,
			#Ok0 = 0.,
			w0 = w0 , 
			wa = wa , 
			Tcmb0 = Tcmb0,
			Neff =Neff, 
			name = name) 

		self._H0 = float(H0)
		self._Ob0   = float (Ob0)
		self._sigmamnu  = float(sigmamnu)
		self._As    = As
		self._sigma8 = sigma8 
		self._ns    = float(ns)

		self._providedtransferfiles = {}
		self._transferredshifts  = None
		if len(TFfilelist ) != len(TFredshiftlist) :
			print "number of files in list must equal number of redshifts in list\n"
			raise  
		if len(TFredshiftlist) > 0:
			for i, z in enumerate(TFredshiftlist) :
				self._providedtransferfiles[z] = TFfilelist[i] 
			self._transferredshifts = sorted (self._providedtransferfiles)
			

	@property 
	def transferredshifts( self ):
		return self._transferredshifts
		
	@property 
	def transferfilesforredshift(self ) :
		return self._providedtransferfiles 
	@property 
	def transferfile ( self , z ) :
		if redshift in transferredshifts :
			return self.transferfilesforredshift( z )  

		else :
			return None

			
	@property 
	def Ob0(self):
		return self._Ob0 

	@property 
	def sigmamnu(self) :
		return self._sigmamnu 


	@property
	def On0 (self) :
			#For now 
		On0 = self._sigmamnu / 94.0 /self.h /self.h
		return On0
	@property
	def Oc0 (self) :
		#print type(self.Om0), self._Ob0 , self._sigmamnu /94.0
		Oc0 = self.Om0 - self._Ob0 - self.On0
		return Oc0


	@property 
	def As(self):
		return self._As

	@property
	def sigma8 (self) :
		return self._sigma8

	@property 
	def ns(self):	
		return self._ns 


	def setamplitude(self, sigma8=None, As=None):

		print sigma8
		print As
		self._sigma8 = sigma8
		self._As = As
	#new methods:
	def summary(self ) :
		h = self.h

		s  = "******************************************************\n"
		s += "======================Cosmology=======================\n" 
		s += " h      = " + str(self.h) +"\n"
		s += "================Relativistic Components===============\n" 
		s += "\Omega_\gamma                  =  " + str(self.Ogamma0)  + "\n" 
		s += "\Omega_nu_0                    =  " + str(self.Onu0)  + "\n" 
		s += "------------------------------------------------------\n"
		relenergydensity  = self.Ogamma0 + self.Onu0 
		s += "sum                            =  " + str(relenergydensity)  + "\n"


		s += "\omega_{nu} (massive)         = "+ str(self.sigmamnu/94.0) + "\n"
		s += "sum of masses (eV )            = "+ str(self.sigmamnu) 
		s += "sum of masses of nu /94.0/h^2  = "+ str(self.sigmamnu/94.0/h/h) + "\n"
		s += "\Omega_{nu} (massive)          = "+ str(self.On0) + "\n"
		s += "\omega_{cdm}                   = "+ str(h*h*self.Oc0) + "\n"
		s += "\omega_{b}                     = "+ str(h*h*self.Ob0) +  "\n"
		s += "\Omega_{cdm}                   = "+ str(self.Oc0) + "\n"
		s += "\Omega_{b}                     = "+ str(self.Ob0) +  "\n"
		s += "------------------------------------------------------\n"
		#s += "--------------------------------\n"
		s += "sum                            =  " + str(self.Oc0 + self.Ob0 + self.On0) +"\n"
		s += "should be the same as \Omega_m =  " + str(self.Om0)  + "\n"
		s += "\Omega_{de}                    =  " + str(self.Ode0) + "\n"
		s += "\Omega_k                       =  " + str(self.Ok0)  + "\n" 
		
		s += "------------------------------------------------------\n"
		s += "sum of components -1.0         =  " + str(self.Oc0 + self.Ob0 + self.On0 +self.Ode0 + relenergydensity -1.) +"\n"
		#s += "=======================\n" 
		s += "=========Primordial Fluctuations =====================\n"
		if self.ns !=None:
			s += "\\n_s                           = "  + str(self.ns) \
				+ "\n" 
		else:
			s += "\\n_s                           = "  + \
				"Not provided \n" 
		if self.sigma8 !=None:
			s += "\sigma_8                       = "                 \
				+ str(self.sigma8) + "\n" 
		else:
			s += "\sigma_8                       = "  + \
				"Not provided \n" 
		if self.As !=None:
			s += "A_s                            = "  + str(self.As) \
				+ "\n" 
		else:
			s += "A_s                            = "  + \
				"Not provided \n" 
		return s

	def growth(self, z ):
		""" 
		returns the linear growth function D0, and the derivative D1 
		of the linear growth function with respect to the scale factor 
		(CHECK) for the cosmology at redshift z

		args:
			z: array-like, make sure they are floats not ints
		returns:
			tuple of D0, and D1 values at the requested redshifts
			
		"""
		import growthfunction as gf

		h = self.h
		omegacb = self.Ob0 + self.Oc0 
		#print "line 93", type(self.On0), type(h)
		omeganuh2 = self.On0 * h * h 
		avals, D, info  = gf.growth(Omegam = omegacb ,
			w0 = self.w0 , wa = self.wa , Tcmb = self.Tcmb0 , 
			h = self.h , Neff = self.Neff , 
			omeganuh2val =  omeganuh2 )


		D , logD = gf.interp_D( z, avals, D[:,0], D[:,1])

		return D, logD 

		
	def growthsourced(self, z ):
		""" 
		returns the linear growth function D0, and the derivative D1 
		of the linear growth function with respect to the scale factor 
		(CHECK) for the cosmology at redshift z

		args:
			z: array-like, make sure they are floats not ints
		returns:
			tuple of D0, and D1 values at the requested redshifts
			
		"""
		import growthfunction as gf

		h = self.h
		omegacb = self.Ob0 + self.Oc0 
		#print "line 93", type(self.On0), type(h)
		omeganuh2 = self.On0 * h * h 
		avals, D, info  = gf.growth(Omegam = omegacb + self.On0,
			w0 = self.w0 , wa = self.wa , Tcmb = self.Tcmb0 , 
			h = self.h , Neff = self.Neff , 
			omeganuh2val =  omeganuh2 )


		D , logD = gf.interp_D( z, avals, D[:,0], D[:,1])

		return D, logD 
if __name__=="__main__":

	f = FCPL ( H0 = 65, Om0 = 0.99 , w0 = -1.0, wa =0.) 
	g = FCPL ( H0 = 65, Om0 = 0.8 , w0 = -1.0, wa =0.) 
	r = FCPL ( H0 = 65, Om0 = 0.3 , w0 = -1.0, wa =0.) 
	f = FCPL ( H0 = 65, Om0 = 0.99 , w0 = -1.0, wa =0., TFfilelist = ["example_data/cmbM000.tf"], TFredshiftlist = [0.0]) 
	#f.setfromHACCinput("example_data/indat.params")
	print "results \n"
	print f.sigmamnu 
	print f.As
	print f.sigma8
	print f.h
	print "H0 ",f.H0
	print f.Om0
	print "get baryon ", f.Ob0
	print "get w0", f.w0
	print "get TCMB ", f.Tcmb0
	print "transfer function file name ", f.transferfilesforredshift[0.0]
	#print "transfer function file name ", f.transferfilesforredshift[0.4]
	#import inspect
	#print inspect.getmembers(f)
	z   = np.arange(0.,5., 0.1)
	plt.plot ( 1.0/(1.0 + z) , f.growth(z)[0] , 'r', label = "Om0 = 0.99")
	plt.plot ( 1.0/(1.0 + z) , g.growth(z)[0] , 'k', label = "Om0 = 0.8")
	plt.plot ( 1.0/(1.0 + z) , r.growth(z)[0] , 'b', label = "Om0 = 0.3")
	plt.plot ( 1.0/(1.0 + z) ,1.0/(1.0 + z) , ls = 'dashed')
	#plt.xlim(0.5,1.0)
	plt.legend(loc ="best")
	plt.show()


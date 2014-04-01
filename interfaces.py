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
		
		w0waCDM.__init__(self, 
			H0 = H0, 
			Om0 = Om0, 
			Ode0 = 1. - Om0,
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



	#new methods:

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


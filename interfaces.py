#!/usr/bin/env python

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
		setvaluesfrom = None , 
		name='FCPL'):
		
		#print "wa ", wa
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
		#self._Oc0   = float (Oc0)
		self._Ob0   = float (Ob0)
		self._sigmamnu  = float(sigmamnu)
		self._As    = As
		self._sigma8 = sigma8 
		self._ns    = float(ns)
		#print self.Tcmb0 
		#print self.w0
		#self._w0 = w0

	@property 
	def Ob0(self):
		return self._Ob0 

	@property 
	def sigmamnu(self) :
		return self._sigmamnu 

	@property 
	def As(self):
		return self._As

	@property 
	def ns(self):	
		return self._ns 

	@property
	def sigma8 (self) :
		return self._sigma8

	@property
	def On0 (self) :
			#For now 
		On0 = self._sigmamnu / 94.0
	@property
	def Oc0 (self) :
		#print type(self.Om0), self._Ob0 , self._sigmamnu /94.0
		Oc0 = self.Om0 - self._Ob0 - self._sigmamnu/94.0 
		return Oc0

	#@property 
	#def w0 (self) :
	#	return self._w0



	def setfromHACCinput(self, inputfile ) :
		from utils import ioutils as io

		haccdict   = io.builddict(inputfile , dictdelim = " " )

		hacc = _hacccosmologydict (haccdict)

		h = hacc['HUBBLE']
		H0 = 100.0*h
		Omega_CDM= hacc['OMEGA_CDM']
		Omega_nu = hacc['OMEGA_NU']
		#print h
		Omega_b = hacc['DEUT']/h/h 
		Omega_m = Omega_CDM + Omega_b + Omega_nu
		self._Om0 =  Omega_m 
		self._Ob0 = Omega_b
		
		self._H0 = H0

		return 0
		
if __name__=="__main__":

	f = FCPL ( H0 = 65, Om0 = 0.3 , w0 = -1.5, wa =0.) 
	#f.setfromHACCinput("example_data/indat.params")
	print f.sigmamnu 
	print f.As
	print f.sigma8
	print f.h
	print "H0 ",f.H0
	print f.Om0
	print "get baryon ", f.Ob0
	print "get w0", f.w0
	print "get TCMB ", f.Tcmb0
	#import inspect
	#print inspect.getmembers(f)
		

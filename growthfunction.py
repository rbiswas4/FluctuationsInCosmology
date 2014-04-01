#!/usr/bin/env python

import sys
import numpy as np

	

	
def omeganuh2a( a , 
	Neff, 
	omeganuh2  = 0. ,
	sigmamass = None , 
	omegagammah2 = 2.471e-5 , 
	method  = "HACC") :

	"""
	returns a tuple of (omeganuh2 (a) , \dot(omeganuh2) (a)) in different 
	approximations specified by the string method. If sigmamass is 
	provided (ie. not None), the mass is assumed  to be equal for each 
	neutrino and overrides the specified values of omeganuh2 using the 
	relation omeganuh2 = sigmamass /94 eV 

	args:
		a	:
		Neff 	:
		omeganuh2 :
		sigmamass: optional, float, defaults to None
			sum of masses of neutrinos assumed in units of eV. 
		method:

	"""
	import utils.typeutils as tu

	if tu.isiterable(a) :
		a  = np.asarray(a) 
	else :
		a = np.asarray([a])

	print a
	sys.exit()
	if sigmamass !=None:
		#mass  = sigmamass /Neff 
		omeganuh2 =  sigmamass / 94.0 


	omeganumasslessh2  = Neff*omegagammah2 
	omeganumasslessh2 *= (7.0/8.0) *(4.0/11.0)**(4.0/3.0)/a**4.0

	omeganumassiveh2a =  omeganuh2 /a**3.0

	if method =="HACC":
		masslessval = -4.0 
		massiveval = -3.0 
		omeganuh2, res = tu.greaterarray( omeganumasslessh2 , omeganumassiveh2a,masslessval, massiveval)
		#print 1./a - 1.0, masslessval , massiveval ,res	, omeganuh2
		dotomeganuh2  = omeganuh2 *res /a 
		
		return omeganuh2, dotomeganuh2 

	else:
		raise ValueError("Method undefined")
		
	 
def HoverH0(z , Omegam , w0 ,wa ,  h , 
	Tcmb = 2.725, Neff= 3.04, omeganuh2val = 0.0, sigmamass = None, terms = False) :
	"""
	returns the value of H/H0 evaluated at a redshift z for a flat 
	cosmology specified by Omegam = Omega_c + Omega_b, a dark energy 
	specified by w0 , wa  and neutrinos

	args:
		z: 
		Omegam : Omega_c + Omega_b at z = 0 
		w0 :
		wa :
		Neff :
		omeganuh2val : float, optional, defaults to 0.0 
			\Omega_{\nu} h^2 for the massive neutrinos. 
		sigmamass: Sum of mass of neutrinos in eV (where all 
			neutrinos are asssumed equal)
			mass of single neutrino = sigmamass/Neff 
			and overrides the omeganuh2val 
		terms : 
	Notes:
		Right now assumes that the masses of neutrinos are all
		equal
	TODO: Think about how to return this later on

	"""
	a = 1.0 /(1.0+ z)
		#Omega values at z = 0 

	Omegagamma = 2.471e-5/h/h * (Tcmb/2.725)**4.0   

		#Bad philosophically, OK approx, but in HACC
	#Omeganumassless = Neff*Omegagamma * (7.0/8.0) *(4.0/11.0)**(4.0/3.0)

	#Calculate the omeganuh2 term at given redshift

	omeganuh2, derivativeh2  = omeganuh2a( a  = a, 
		Neff = Neff, 
		omeganuh2 = omeganuh2val , 
		sigmamass = sigmamass, 
		omegagammah2 = Omegagamma * h *h ,
		method =  "HACC")
	
	omeganuh20 , derivative0 = omeganuh2a( a  = 1.0, 
		Neff = Neff, 
		omeganuh2 = omeganuh2val , 
		sigmamass = sigmamass, 
		omegagammah2 = Omegagamma * h *h ,
		method =  "HACC")


	Omegade    = 1.0 - Omegam  - Omegagamma - omeganuh20/h/h
		#Relative density at different redshift
	rhode = Omegade * (1.0 + z )**(3.0 *(w0 + wa + 1.0))
	rhode *= np.exp(-3*wa*z/(1. + z))

	rhom  =  Omegam * (1.0 + z ) ** 3.0 

	rhogamma = Omegagamma * (1.0 + z )**4.0

	rhonu = omeganuh2 / h / h 
	#print rhonu

	#rhonumassless = Omeganumassless * (1.0 + z )**4.0
	
		#HoverHzsq
	Esq =  rhom + rhode + rhogamma + rhonu

		#Calculate derivatives with respect to a
		#from the stress tensor conservation of 
		#individual components

	drhom  = -3.0 * (1.0+ 0.) *rhom / a 
	drhode = - 3.0 *(1.0 + w0 + wa *(1.0- a)) *rhode/ a 
	drhogamma = -3.0 *(1.0+ 1.0/3.0)*  rhogamma /a 
		#Separate: depends on approximation as w is not known
	drhonu  = derivativeh2/h/h 

	if terms :
		return np.sqrt(Esq) , drhom , drhode  , drhogamma , drhonu
	else:
		return np.sqrt(Esq)
	

def dlnHda (a , Omegam , w0, wa, h, Tcmb = 2.725, Neff = 3.04, omeganuh2val = 0.0) :

	
	z = 1./a - 1.
	H, drhom , drhode  , drhogamma , drhonu = HoverH0 ( z , 
		Omegam , w0, wa ,  
		h, Tcmb = Tcmb, omeganuh2val = omeganuh2val , 
		Neff = Neff, terms = True) 

#	drhom  = -3.0 * rhom / a 
#	drhode = - 3.0 *(w0 + wa + 1.0 - a*wa) *rhode/ a 
#	drhogamma = -4.0 *  rhogamma /a 
#	drhonu  = derivative

	#drhonumassless = -4.0 *  rhonumassless /a 
	dH2da = (drhom + drhode + drhogamma + drhonu )
	
	return dH2da/H/H /2.0
	

def Ddelta ( Delta , a , Omegam, w0, wa , h , Tcmb = 2.725 , Neff = 3.04, omeganuh2val = 0.0, terms = False) :

	z = 1.0/a -1. 

	HoverH0sq = HoverH0(z , Omegam ,  w0, wa , h = h, Tcmb = Tcmb, Neff = Neff, omeganuh2val = omeganuh2val) 
	HoverH0sq = HoverH0sq* HoverH0sq

	#print "in Ddelta" 
	#print a , Omegam , w0 , wa , h , Tcmb , Neff
	dlnHdaval = dlnHda (a , Omegam , w0, wa , h = h , Tcmb = Tcmb, Neff = Neff, omeganuh2val = omeganuh2val)  
	
	DY0 = Delta[1] 

	DY1coeffDelta1 = -( 3.0/ a + dlnHdaval ) 
	#print DY1
	#DY1 *= Delta[1]  						
	#print DY1
	#print Omegam , Delta[0], HoverH0sq
	DY1coeffDelta0  = 3.0*Omegam /2.0 / a**5.0 / HoverH0sq 
	DY1 = DY1coeffDelta0 * Delta[0] + DY1coeffDelta1*Delta[1] 
	#print DY1

		#The useful print statement
	#print a , DY1coeffDelta0 , DY1coeffDelta1 ,DY0 , DY1
	if terms:
		return [DY0, DY1] , DY1coeffDelta1 , DY1coeffDelta0
	else:
		return [DY0 , DY1 ] 

def growth ( Omegam , w0 , wa , Tcmb , h  , Neff , omeganuh2val ):

	from scipy.integrate import odeint

	ainit = 1.0e-7
	afinal = 1.0 
	steps = 100000

		###Calculate Initial Delta, so that 2nd derivative stays 
		###manageable at 0 derivative at high redshifts
	zinit = 1.0/ainit - 1.0 
	H0overH = 1.0/HoverH0(zinit, Omegam , w0 , wa , h ,Tcmb , Neff , omeganuh2val = omeganuh2val)
		# Order of magnitude of source term	
	orderofmag = H0overH*H0overH *ainit**(-5.0)
	logainit= np.log(ainit)
	logafinal = np.log(afinal)
		#Initial conditions. Note that the value of Delta is 
		#artificial. Set to inverse of order of magnitude of source
		#term. 
		#We will normalize later so that D(a=1) = 1.
	#Dinit = np.array([1.0/orderofmag , 0.0 ])
	Dinit = np.array([ainit , 0.0 ])
	#print Dinit
	#avals = np.logspace(logainit , logafinal ,steps)	
	avals = np.linspace(ainit, afinal, num = steps)
		#Not normalized	
	growthun , info = odeint(Ddelta, Dinit, avals, full_output =1, mxstep = 100000,args=(Omegam, w0, wa, h, Tcmb , Neff, omeganuh2val)) 

		# Normalization
	growth = growthun
	growth[:,0] = growthun[:,0]/growthun[-1,0]
	growth[:,1] = growthun[:,1]/growthun[-1,0]

	
	return avals, growth ,info 

def interp_D(zseek, ascale, D0, D1):
	"""
	interpolate 1D, using numpy, the given values of the growth 
	function D0 and its derivative D1 at given values of the scale 
	factor ascale to return a tuple of the growth function Dseek, 
	and logDseek at the requested values of redshift zseek.  

	args:
		zseek  :

		ascale :

		D0     :
		
		D1     :

	returns:
		tuple (Dseek , logDseek ) , where Dseek is the interpolated 
		values of D0 at zseek, and logDseek, is the interpolated
		values of D1 at zseek.

	status:
		Written by Suman Bhattacharya, Dec 2013
	"""

	aseek    = 1/(1+zseek)
	Dseek    = np.interp(aseek, ascale, D0, left = np.nan , right = np.nan)
	logDseek = np.interp(aseek, ascale, D1, left = np.nan , right = np.nan)

	return (Dseek, logDseek)

if __name__=="__main__":

	from astropy.cosmology import set_current
	from astropy.cosmology import get_current
	from astropy.cosmology import w0waCDM  

	import numpy as np
	#import matplotlib.pyplot as plt



	#m = w0waCDM(H0=71.0, Om0 = 0.25, Ode0 = 0.75, w0 = -0.80 , wa = 0.30 , Tcmb0 = 0., Neff =0.) 
	#mn = w0waCDM(H0=71.0, Om0 = 0.25, Ode0 = 0.75, w0 = -0.80 , wa = 0.30 , Tcmb0 = 2.725, Neff =3.04) 

	#set_current(m)
	avals, D , info  = growth(Omegam =0.25, w0 =- 1.00, wa =0.00, Tcmb =2.725, h =0.71, Neff  = 3.04,omeganuh2val = 0.01 ) 
	#avals, Dn , infon  = growth(Omegam =0.25, w0 =- 1.0, wa =0.0, Tcmb =2.725, h =0.71, Neff  = 3.04, omeganuh2val = 0.0) 
	#avals, D , info  = growth(Omegam =0.3, w0 =- 1.00, wa =0.00, Tcmb =0, h =0.71, Neff  = 0.00, ) 
	zseek= 1.0
	zseek = np.arange(0.0, 2.0, 0.1)
        Dseek, logDseek=interp_D(zseek, avals, D[:,0],D[:,1])
        print np.array([zseek, Dseek, logDseek]).T

#!/usr/bin/env python


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

		BUGFIX: unittype ought to be kgperm3 in if loop, but was 
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


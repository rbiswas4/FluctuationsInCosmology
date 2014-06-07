#!/usr/bin/env python
import numpy as np

def __SP (
	A0 , 
	alpha  ,  
	z = 0.) :
	"""

	"""


	A = A0 / (1. + z )**alpha
	return A 

def fsigmaTinker( sigma, 
   A = 0.186 , 
   a = 1.47 , 
   b = 2.57 , 
   c = 1.19): 

	x =  sigma
 	return A*(((x/b)**(-a))+1)*np.exp(-c/x**2)
def fsigmaMICE(
	sigma , z ) :
	"""
	"""
	z = np.asarray(z)

	A = 0.58 * (1.0 + z )**(-0.13)
	a = 1.37 * (1.0 + z ) ** (-0.15)
	b = 0.3  * (1.0 + z )**(-0.084)
	c = 1.036 * (1.0 + z) **(-0.024)

	if z > 1.0 :
		return np.ones(len(sigma))* np.nan

 	return A * ( sigma**(-a) + b ) *np.exp( - c/sigma/sigma )
#def fsigmaBhattacharya(
#	Mass , 
#	deltac = 1.674 , 
#	z = 0.0 ,
#	A0  = 0.333 ,
#	a0 = 0.788 , 
#	p0 = 0.807 , 
#	q0 = 1.795 , 
#	alpha1 = 0.11 ,  
#	alpha2 = 0.01 , 
#	alpha3 = 0.0 , 
#	alpha4 =  0.0,  
#	Mlow = 6e11 , 
#	Mhigh = 3e15) :
#
#	"""
#	"Produces the fit from Bhattacharya etal. in Table 4 of 
#	http://arxiv.org/pdf/1005.2239.pdf
#	f(sigma ) =  A(z) (1.0 +(\sigma^2/\delta_c^2/a)^p)\frac{2}{\pi}\exp( -\frac{a \delta^2_c}{2\sigma^2}) 
#	
#	
#	args:
#
#	returns: 
#	"""	
#	import psutils as psu
#
#	#return only values within mass range where fit is advertised
#	Ml = Mass > Mlow  
#	Mh = Mass < Mhigh   
#
#	Mass = Mass[Ml & Mh ]
#
#	sigmaM = psu.sigmaM(M= Masses, ps = ps, cosmo = cosmo, z ) 
#
#	return __fsigmaBhattacharya ( sigmaM , 
#		deltac = 1.674 ,
#	A0 = 0.333 , 
#	a0 = 0.788 , 
#	p0 = 0.807 , 
#	q0 = 1.795 , 
#	alpha1 = 0.11 ,  
#	alpha2 = 0.01 , 
#	alpha3 = 0.0 , 
#	alpha4 =  0.0,  
#	Mlow = 6e11 , 
#	Mhigh = 3e15)
def __fsigmaBhattacharya (  
	sigma ,
	deltac = 1.686 , 
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
	Mhigh = 3e15) :

	"""Produces the fit from Bhattacharya etal. in Table 4 of 
	http://arxiv.org/pdf/1005.2239.pdf


	f(sigma ) =  A(z) (1.0 +(\sigma^2/\delta_c^2/a)^p)\frac{2}{\pi}\exp( -\frac{a \delta^2_c}{2\sigma^2}) 
	"""


	#Check mass range
	deltac2 = deltac*deltac
	sigma2 = sigma*sigma

	A = __SP (A0 , alpha1, z )
	a = __SP (a0 , alpha2,  z )
	p = __SP (p0, alpha3, z )
	q = __SP (q0, alpha4, z )  

	#print A , a, p, q 
	

	f  = A * np.sqrt(2.0 / np.pi) 
	#print f
	xx = ( 1. + ( sigma2/a / deltac2 )**p )
	#print xx
	f *= xx
	yy = (np.sqrt(a * deltac2 /sigma2 ))**q
	#print yy 
	f *= yy  
	zz = np.exp ( - a * deltac2/sigma2 /2.0 ) 
	#print zz 
	f *= zz 

	return f

	

def __fsigmaShethTormen ( 
	sigma , 
	deltac  = 1.686 ,
	z = 0.0 ,
	A0 = 0.322 , 
	a0 = 0.75 , 
	p0 = 0.3 , 
	q0 = 1.0 , #Defines Sheth Tormen  
	alpha1 = 0.0 ,  
	alpha2 = 0.0 , 
	alpha3 = 0.0 , 
	alpha4 =  0.0,  
	Mlow = 6e11 , 
	Mhigh = 3e15) :
	
	

	return __fsigmaBhattacharya( sigma , deltac , z , A0 , a0 , p0 , q0 , 
		alpha1, alpha2 , alpha3 , alpha4 , Mlow , Mhigh )


def bias_Bhattacharya( 
	sigma , 
	deltac  = 1.674 ,
	z = 0.0 ,
	a0 = 0.75 , 
	p0 = 0.3 , 
	alpha2 = 0.01 , 
	alpha3 = 0.0 , 
	Mlow = 6e11 , 
	Mhigh = 3e15) :
	"""

	"""


	nu = z

	
	bias =  1.0 + (a*nu - 1.0)/ deltac 
	bias += 2.0*p / deltac /(1.0 + (a*nu)**p)

	return bias
if __name__ == "__main__":


	import numpy as np 
	import matplotlib.pyplot as plt 



	oneoversigma = np.arange(0.1, 3, 0.01) 

	sigma = 1.0 / oneoversigma  
	fsigmaB = __fsigmaBhattacharya (sigma )
	fsigmaST = __fsigmaShethTormen (sigma )
	plt.loglog ( oneoversigma , fsigmaB, label = "Bhattacharya etal.") 
	plt.loglog ( oneoversigma , fsigmaST, label = "Sheth-Tormen") 
	plt.xlabel(r'$1/\sigma$')
	plt.ylabel(r'$f(\sigma)$') 
	plt.legend(loc = "best")
	plt.tight_layout()
	plt.figure()
	plt.plot(oneoversigma , fsigmaST/fsigmaB, label= "Sheth-Tormen")
	plt.xlabel(r'$1/\sigma$')
	plt.ylabel(r'$f(\sigma)/f_{Bh}(\sigma)$')
	plt.xlim(0.4, 2.8)
	plt.ylim(0.7,1.2)
	plt.legend(loc="best")
	plt.tight_layout()
	plt.show()


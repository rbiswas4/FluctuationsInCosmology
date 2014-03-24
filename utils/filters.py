#!/usr/bin/env python
# Set of isotropic filters to use in calculations


import numpy as np
import typeutils as tu

def Wtophatkspace(kR, R = None ) :
	"""Returns the Fourier Transform of the real space top hat window 
	function which is 1 when r < R, and 0 otherwise.
	args:
		kR: unit less quantity
	returns: 
		array like, window function in k space
	"""
	filt  = 3*(np.sin(kR) - (kR)*np.cos(kR))/(kR)**3
	return filt


def Wtophatkspacesq (kR , R = None):

	filt  = Wtophatkspace(kR, R )

	return filt*filt
def dWtophatkspacedR ( kR  ,R ) :


	"""
	Returns the derivative of the Fourier Transform of the real space 
	top hat window function which is 1 when r < R, and 0 otherwise with 
	respect to R. 

	args:
		kR : unit less quantity. 
		R : Radius 
	returns:
		array like derivative of Wtophatkspace with respect to R

	"""
	filt = 3.0 *(  np.sin(kR) / kR  - Wtophatkspace (kR )) 
	if tu.isiterable(R):

		x = filt.transpose()/R
		#returns an array whose elesements are k values for each R
		return x.transpose()

		
	else:

		return  filt /R 

def dWtophatkspacedlnR ( kR  ) :


	"""
	Returns the derivative of the Fourier Transform of the real space 
	top hat window function which is 1 when r < R, and 0 otherwise with 
	respect to R. 

	args:
		kR : unit less quantity. 
	returns:
		array like derivative of Wtophatkspace with respect to R

	"""
	filt = 3.0 *(  np.sin(kR) / kR  - Wtophatkspace (kR )) 
	#if tu.isiterable(R):
	#	x = filt.transpose()/R
	#	#returns an array whose elesements are k values for each R
	#	return x.transpose()

		
	#else:

	#	return  filt /R 
	return filt 

def dWtophatkspacesqdR( kR , R ) :


	w1 = dWtophatkspacedR (kR, R) 
	w2 = Wtophatkspace(kR , R)
	#print "w1,w2", w1 , w2
	return 2.0*w1 * w2 
def dWtophatkspacesqdlnR( kR , R = None) :


	w1 = dWtophatkspacedlnR (kR) 
	w2 = Wtophatkspace(kR )
	#print "w1,w2", w1 , w2
	return 2.0*w1 * w2 


if __name__ == "__main__":

	import matplotlib.pyplot as plt
	import numpy as np

	k = np.logspace(-6,5,100)
	#R = np.arange(1,8,1.)
	#kr  = np.outer (k,R) 

	plt.plot(k , Wtophatkspace(k*2.0)**2.0, label = "R= 2.0")
	plt.plot(k , Wtophatkspace(k*5.0)**2.0, label = "R =5.0")
	plt.plot(k , Wtophatkspace(k*8.0)**2.0, label = "R= 8.0")
	plt.legend(loc= "best")
	plt.ylabel("Wtophatkspace(k*R)")
	plt.xlabel(r'$k h^{-1} Mpc$')
	plt.xscale('log')
	plt.show()

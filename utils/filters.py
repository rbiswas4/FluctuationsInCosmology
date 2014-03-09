#!/usr/bin/env python
# Set of isotropic filters to use in calculations


import numpy as np

def Wtophatkspace(kR) :
	"""Returns the Fourier Transform of the real space top hat window 
	function which is 1 when r < R, and 0 otherwise.
	args:
		kR: unit less quantity
	returns: 
		array like, window function in k space
	"""
	return 3*(np.sin(kR) - (kR)*np.cos(kR))/(kR)**3


def dWtophatkspacedR ( k, R  ) :


	"""
	Returns the derivative of the Fourier Transform of the real space 
	top hat window function which is 1 when r < R, and 0 otherwise with 
	respect to R. 

	args:
		k : unit less quantity. 
		R : Radius 
	returns:
		array like derivative of Wtophatkspace with respect to R

	"""

	return  3.0 (  Sin(k*R) / R / k  - Wtophatkspace (k*R )) /R 


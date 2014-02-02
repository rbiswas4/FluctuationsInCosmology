#!/usr/bin/env python

import collections
import numpy as np
def greaterarray( a1 , a2 ,a1val =None, a2val = None) :

	"""
	Given two arrays a1 and a2 of equal length, return an array with the 
	ith element of a1 if a1[i] > a2[i] and a2[i] otherwise. 

	"""
	if a1val ==None:
		a1val = 0 

	if a2val == None:
		a2val = -1

	if not isiterable (a1) :

		if isiterable(a2) :

			raise ValueError("a2 is array, and a1 is scalar")

		if a1 > a2 :
			return a1, 0 
		else:
			return a2, -1 

	else:
		if not isiterable(a2) :
			raise ValueError("a1 is array, and a2 is scalar")
	
		if len(a2) == len(a1) :
			res = np.zeros(len(a1))
			for i in range(len(a1)):
				if a1[i] > a2[i] :
					res[i] = a1val 
				else :
					res [i] = a2val 
			return np.fmax(a1,a2), res  

def findtype(s, makeintfloats = False):

	"""
	Return a tuple with the data type of the string along with 
	the string converted to the data type.

	args:
		s : mandatory, string
	returns:
		tuple: (type, sprime)
			type is either int, float or str
			sprime is the quantaty s converted into its type.
	example usage:
		t, s = findtype(s) 
	status:
		seems to be working,
		R. Biswas, Sun Mar 24 21:40:53 CDT 2013
		Copied from ioutilst as is
		R. Biswas, Tue Oct 22 14:47:48 CDT 2013
		 
	"""
	try:
		int(s)
		if makeintfloats:
			return 'f4', float(s)
		else:
			return 'i8' , int(s) 
	except ValueError:
		pass
	try:
		float(s)
		return 'f4' , float(s) 
	except ValueError:
		pass
 
	return "a20", s
def isiterable ( var ) :
	"""
	returns True if var is an iterable quantity, False otherwise
	
	args:
		var: mandatory, anytype
	returns : bool 
		True if iterable, False otherwise
	example_usage:
		>>> if io.isiterable(k) :	
		>>>	k = np.asarray(k)

	status:
		tested, and works. 
		R. Biswas, Wed Oct  9 11:05:12 CDT 2013
	"""
	return isinstance(var , collections.Iterable)

def hasanyattrs (object ,returnlist = False):

	"""
	Checks if the object has any attributes (rather than a particular 
	attribute as implemented in python 2.7.5 function hasattr ) 

	args: 
		object : name of object 

	returns: 
		tuple, ( Bool, listof attributes) 
	example usage:
		>>> import typeutils as tu 
		>>> print tu.hasanyattrs(o ,True)
		
	status: 
		Tested and found to work as expected, 
		R. Biswas, Sun Oct 13 17:01:45 CDT 2013
		
	"""

	hasanyattrs = False 
	lst = object.attrs.items()
	if len(lst) !=0:
		hasanyattrs = True 

	return (hasanyattrs , lst )

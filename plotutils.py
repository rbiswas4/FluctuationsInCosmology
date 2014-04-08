#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def settwopanel( height_ratios = [1.0,0.3], 
	width_ratios = [1.,0.] ,
	padding = None, 
	setoffset = None, 
	figsize  = None ):

	"""
	returns a figure and axes for a main panel and a lower panel for 
	showing differential information of overplotted quantities in
	the top panel. 
	args	:
		height_ratios: list of floats, optional defaults to 
				[1.0, 0.3]
			height ratio between the upper and lower panel 

		width_ratios :list of floats, optional defaults to 
				[1.0, 0.0]
			width ratio between the left and right  panel 
		
		figsize: 

	returns :
		figure object , ax0 (axes for top panel) , and ax1 
			(axes for lower panel)

	usage   :
		>>> myfig,  myax0 , myax1 = settwopanel ( )
		>>> myax0.plot( x,  y) 
		>>> myax1.plot(x, x)
		>>> myfig.tight_layout()

	status  :

		tested by 
		R. Biswas, Fri Feb 21 00:52:55 CST 2014
	"""
	import matplotlib.ticker as ticker
	majorformatter = ticker.ScalarFormatter(useOffset =False)

	if figsize == None:
		fig = plt.figure()
	else:
		fig = plt.figure(figsize = figsize)
	

	gs = gridspec.GridSpec ( 2, 1, width_ratios = width_ratios , height_ratios = height_ratios)

	ax0 = plt.subplot(gs[0]) 
	ax1 = plt.subplot(gs[1])
	ax0.set_xticklabels("",visible = False)
	ax1.yaxis.set_major_formatter(majorformatter)

	hpad  = 0.0 
	#gridspec.update(hpad = hpad)
	return fig , ax0 , ax1


if __name__ == "__main__":


	x = np.arange(0,10,0.1)
	y = x * x

	
	myfig,  myax0 , myax1 = settwopanel ( )

	myax0.plot( x,  y) 
	myax1.plot(x, x)
	
	myfig.tight_layout()
	plt.show()

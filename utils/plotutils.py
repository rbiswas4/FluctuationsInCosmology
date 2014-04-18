#!/usr/bin/env python

import typeutils  as tu
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def drawxband(refval , 
	xlims = None, 
	color = 'gray',
	bandwidths = [-0.1,0.1], 
	ua = None ):
	
	"""
	draw a shaded band of color color of certain width around a single 
	reference value refval
	args:
		refval    :	float , mandatory
			scalar (single) reference value, 
			Example: refval  = 1.0
			
		xlims     : list of two floats, optional, defaults to None
			min and max x values through which the band will 
			be drawn
			if None, these values are set from the limits of the 
			supplied axessubplot. If axessubplots is not supplied
			this will raise an error
				
		color     : python color style, optional,defaults to 'gray' 
			color of the band
		bandwidths : list of two floats, optional default to [-0.1, 0.1]
			width of the band to be shaded. 
		ua	  : optional, axes.subplot object on which to draw, 
			defaults to None
			axes object to use
	returns:
		axes object
	example usage:
		>>> #To use an existing axes subplot object
		>>> drawxband (refval = 60., bandwidths = [-20.,20.], color = 'green', ua = myax0)	
		>>> #No axes object is available
		>>> drawxband (refval = 60., xlims = [4., 8.] , bandwidths = [-20.,20.],color = 'green')	

	status: Tests in plot_utils main. 
		tested R. Biswas, Fri Apr 18 08:29:54 CDT 2014 

	"""

	if ua == None:  
		ua = plt.gca()
	else:
		plt.sca(ua)
		xl ,xh = ua.get_xlim()
		if xlims == None:
			xlims = [ xl , xh ] 


	#really started this conditional statement to make it general, 
	#but have not finished, will always go to the else condition
	#in current implementation
	if xlims == None:
		if not tu.isiterable(refval):
			raise ValueError("To draw band supply either refval as numpy array or xlims over which it should be drawn")
				
	else:
		xvals = np.linspace(xlims[0], xlims[1],2)

			#refvals = refval *np.ones(len(xvals) 



	#plot reference value
	ua.axhline(refval,color ='k',lw =2.0)
	#draw band
	ua.fill_between (xvals , refval + bandwidths[0], refval + bandwidths[1],
		color  = color , alpha =  0.25) 

	return ua
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
	drawxband (refval = 60.,bandwidths = [-20.,20.],color = 'green', ua = myax0)	
	#myax1.plot(x, x)
	
	myfig.tight_layout()

	plt.figure()
	plt.plot(x, y)
	drawxband (refval = 60., xlims = [4., 8.] , bandwidths = [-20.,20.],color = 'green')	
	
	
	plt.show()

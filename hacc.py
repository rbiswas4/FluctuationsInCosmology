#!/usr/bin/env python
#Requirements:
#cosmoutils 
#
from utils import ioutils as io
import numpy as np
import matplotlib.pyplot as plt
import sys

#example input file for HACC
haccindat      =  "example_data/indat.params"
#haccdictfile   = io.builddict(haccindat , dictdelim = " " )
def haccredshiftsforstep(stepnumber ,
	zin = 200. , 
	zfinal  = 0. ,
	numsteps = 500) :

	ain  =  1.0/(1.0 + zin) 
	afinal  =  1.0/(1.0 + zfinal) 
	astep =  (afinal - ain)/numsteps 

	aatstep = astep * stepnumber + ain 
	z = 1.0/aatstep - 1.0 

	return z 

def particlemass(
	numparticlescuberoot ,
	lengthcuberoot, # cube root of vol in h^{-1} Mpc
	Omega , 
	h = 1.0) :

	"""returns the mass assigned to N-Body particles of the particular 
	species given the number of particles in the simulations, the length 
	scale of simulation, and the relative density of the species wrt 
	the critical density and the hubble parameter h = H0/100.

	args:
		numparticlescuberoot : float, mandatory
			cube root of the number of particles representing the
			species in the simulation, 
			ie. numparticles= numparticlescuberoot^3 

		lengthcuberoot : float, mandatory 
			cube root of the length scale of the (assumed to be 
			cubic) simulation in units of h^{-1} Mpc. Thus the 
			vol of the simulation is (lengthcuberoot*h)^3 Mpc^3

		Omega : float, mandatory
			relative density of the species in units of critical 
			density at z =0 . ie rho_{species} /rho_{crit}(z =0). 

		h     : float, optional, defaults to 1.0 
			H0/h used in simulation
	returns:
		mass assigned to particles of that species in the N-Body
		simulation in units of solar masses. 		

	example usage:
		>>> print particlemass (numparticlescuberoot = 256 , 
			lengthcuberoot = 512 ,
			Omega = 0.25 , 
			h = 0.7)
		>>> 93265221934.8

	status:
		Not tested, 
	todo:
		Test with one known case in HACC

	"""
	from astropy import units as u
	from cosmoutils import critdensity


	#actual rhocrit0 = rhocrit0 * h^2
	rhocrit0 = critdensity(h = 1., unittype = 'solarmassperMpc3') 
	#units of solar masses/ Mpc^3  / h^2
	#print "RHOCRIT", rhocrit0
	rho = Omega * rhocrit0 

	numparticles = numparticlescuberoot**3.0 
		#Get num density in m^{-3}
	#vol = (lengthcuberoot * h) **3.0  * (u.Mpc.to(u.m))**3.0
	#numdensity =  numparticles / vol 

	#number density in number/ (h^{-1} Mpc) ^3
	numdensity = (numparticlescuberoot/lengthcuberoot)**3.0 

		#particle mass in h^{-1} solar masses 
	particlemass = rho / numdensity 
	
		#particle mass in solar masses
	#convtosolmass =  u.kg.to(u.solMass)	
	return particlemass #*convtosolmass 


def simdict( haccinputdict) :

	simvariables = ['RL' , 'NP' , 'Z_IN', 'N_STEPS', 'Z_FIN']
	haccsim = {}
	for key in haccinputdict :
		if key in simvariables :
			#print key, float(haccinputdict[key])
			haccsim[key] = float( haccinputdict[key]) 

	return haccsim 

def hacccosmologydict(haccinputdict ) :

	cosmologyvariables =  ['HUBBLE', 'OMEGA_CDM','DEUT','OMEGA_NU' ,'W_DE','WA_DE', 'SS8', 'NS']
	hacccosmologydict = {}
	for key in haccinputdict :
		if key in cosmologyvariables :
			hacccosmologydict[key] = float(haccinputdict[key])

	return hacccosmologydict 

def _hacccosmo(inputfile ) :
		from utils import ioutils as io

		haccdict   = io.builddict(inputfile , dictdelim = " " )

		#print haccdict
		#print haccdict['OMEGA_NU']
		#sys.exit()
		hacc = hacccosmologydict (haccdict)

		h = hacc['HUBBLE']
		H0 = 100.0*h

		Omega_CDM= hacc['OMEGA_CDM']
		Omega_nu = hacc['OMEGA_NU']
		Omega_b = hacc['DEUT']/h/h 
		Omega_m = Omega_CDM + Omega_b + Omega_nu

		w0 = hacc['W_DE']
		wa = hacc['WA_DE']
		sigma8  = hacc['SS8']
		ns = hacc['NS']

		return (H0 , Omega_m , Omega_nu, Omega_b , w0, wa , sigma8, ns) #check 
		


 
class haccsim (object ) :

	def __init__(self, indatfile , 
		name  = None , 
		TKFile = None ) :



		if name == None:
			self.name = ""
		else:
			self.name = name 
		haccdict = io.builddict(indatfile , dictdelim = " " )

		#set cosmo:
		from interfaces import FCPL 
		H0 , Omega_m , Omega_nu , Omega_b, w0, wa, sigma8  , ns =  _hacccosmo(indatfile) 
		sigma_mnu = Omega_nu*94.0 
		#print sigma_mnu 
		#sys.exit()
		
		self.cosmo = FCPL ( H0 = H0 , Om0 = Omega_m , Ob0 = Omega_b , 
			sigmamnu = sigma_mnu , w0 = w0, wa = wa, 
			sigma8 = sigma8  , ns = ns)


		#set sim properties
		haccsimdict = simdict (haccdict )

			#set massresolution
		Omegacb0 = self.cosmo.Ob0 + self.cosmo.Oc0
		self._boxsize = haccsimdict ['RL']
		self._numparticlescuberoot = haccsimdict['NP']
		self._particlemass = particlemass(Omega = Omegacb0 , 
			lengthcuberoot = self._boxsize,
			numparticlescuberoot = self._numparticlescuberoot ,
			h = self.cosmo.h ) 

			#set time stepping info
		self._zin = haccsimdict["Z_IN"]
		self._zf  = haccsimdict["Z_FIN"]
		self._numtimesteps = haccsimdict["N_STEPS"]

	@property 
	def numparticles (self) :
		npcuberoot = self._numparticlescuberoot 
		return npcuberoot * npcuberoot * npcuberoot 
		
	@property
	def particlemass (self) :
		return self._particlemass 

	@property 
	def boxsize (self) :
		return self._boxsize 

	@property 
	def simvolume (self) :
		l =  self._boxsize
		return l*l*l

	@property 
	def initialredshift (self ) :
		return self._zin 

	@property 
	def numtimesteps (self ) :
		return self._numtimesteps 

	@property 
	def finalredshift (self ) :
		return self._zf 

	#Methods 

	def steptoredshift(self , step ) :

		zin =  self.initialredshift  
		#print zin
		zf  =  self.finalredshift 
		#print zf
		numsteps = self.numtimesteps 

		return  haccredshiftsforstep(step ,
			zin = zin , 
			zfinal  = zf ,
			numsteps = numsteps)
	
	def summary (self , logfilename = None) :


		s  = "SUMMARY OF simulation, name = "+ self.name  
		s += "\n\n========================================\n\n"
		s += "Cosmology\n====================\n" 
		s += "massive neutrino energy density " + str(self.cosmo.sigmamnu/94.0) + "\n"
		s += "massive neutrino energy density " + str(self.cosmo.On0) + "\n"
		s += "Simulation Properties\n=====================\n\n"

		s += "Simulation Volume ("+ "{:.2e}".format(self._boxsize) +r'$)^3 h^{-3} Mpc^3$'+" \n"  
		s += "Particle Mass = "+ "{:.2e}".format(self.particlemass)+"r'$h^{-1} M_\odot'$\n" 
		if logfilename != None:

			f = open(logfilename, "w") 
			f.write(s) 
			f.close()

		return s
		

def tests() :

	if len(sys.argv) <2:
		haccindat = "example_data/indat.params"
		#"/media/star1/simulations/Neutrinos/M000n1/output"
	else:
		haccindat = sys.argv[1]
		
	hacc = haccsim ( haccindat, name = "M000")

	print hacc.summary()
	#print hacccosmologydict(haccindat)
	sys.exit()
	print hacc.cosmo.Om0
	print "particle mass in sim is ", hacc.particlemass
	print "mass res in sim is ", "{:.2e}".format(hacc.particlemass )
	print type(hacc.initialredshift)
	print hacc.finalredshift
	#print hacc.steptoredshift ( 499) 
	print hacc.steptoredshift ( 161) 


	print hacc.cosmo.growth(1.0)
	z  = hacc.steptoredshift ( np.arange(1,500,1))
	plt.plot ( 1.0/(1.0 + z) , hacc.cosmo.growth(z)[0] , 'r')
	plt.plot ( 1.0/(1.0 + z) ,1.0/(1.0 + z) , ls = 'dashed')
	plt.xlim(0.5,1.0)
	plt.show()


if __name__=="__main__":

	tests()


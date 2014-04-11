
# In[1]:

import sys
sys.path.append("/Users/rbiswas/src/interfacecosmology")


import halomass as hm
from interfaces import FCPL
import psutils as psu
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import camb_utils.cambio as cio
import psutils as psu
from interfaces import FCPL
import massfunctions as mf
import utils.plotutils as pu
import halomass as hmf
import sys
import matplotlib.gridspec as gridspec


# In[2]:

# INPUT FILES (CAMB outputs)
dirn  = "./"
#M000 (LCDM)
tkfile  = dirn +"example_data/M000/m000n0_transfer_fin_out.dat"
pkfile = dirn + "example_data/M000/m000n0_matterpower_fin.dat"
#M000n1 (same LCDM, with some fraction of CDM replaced by massive neutrinos)
ntkfile = dirn +"example_data/M000n1/m000n1_transfer_fin_out.dat"
npkfile = dirn +"example_data/M000n1/m000n1_matterpower_fin.dat"


# In[3]:

#Cosmological Models:
#               Same cosmology represented differently
M000 = FCPL(H0 = 71. ,Om0  = 0.26479, Ob0 = 0.044793, ns = 0.963, As = 2.16e-9)
M000s  = FCPL(H0 = 71., Om0 = 0.26479, Ob0 = 0.044793, ns = 0.963, sigma8 = 0.8)
#
M000n1 = FCPL(H0 = 71., Om0 = 0.26479, Ob0 = 0.044793, ns = 0.963, sigma8 = 0.8, sigmamnu=0.94) 


# In[4]:

#powerspectrum 
#use power spectrum from CAMB power spectrum output
psfrompk = psu.powerspectrum(koverh = None, asciifile = pkfile)
#use matter power spectrum from CAMB transfer function output, with As
psfromtkm = psu.powerspectrum(koverh = None, pstype ="matter", asciifile =tkfile, cosmo = M000)
#use matter power spectrum from CAMB transfer function output, with sigma8
psfromtkms = psu.powerspectrum(koverh = None, pstype ="matter", sigma8type = "matter", asciifile =tkfile, cosmo = M000s)
#use cb power spectrum from CAMB transfer function output, normalized to matter
psfromtkcbs = psu.powerspectrum(koverh = None, pstype ="cb", sigma8type = "matter", asciifile =tkfile, cosmo = M000s)
#use cb power spectrum from CAMB transfer function output, normalized to cb
psfromtkcbscb = psu.powerspectrum(koverh = None, pstype ="cb", sigma8type = "cb", asciifile =tkfile, cosmo = M000s)
#
#psfrompkn = psu.powerspectrum(koverh = None, asciifile  = npkfile)
#psfromtknm = psu.powerspectrum(koverh = None, asciifile  = ntkfile, pstype ="matter", cosmo = M000n1)


# In[5]:

simdir = "/Users/rbiswas/doc/projects/NeutrinoMF/doc/files/paperfiles/"
fname0 = simdir + "m000.499.fof_mf"


# In[13]:

#M000_z0 = np.loadtxt(simdir + "m000.499.fof_mf")


# In[15]:

#MinhinvMsun = M000_z0[:,0]


# In[25]:

def fracPoissonerrors(num ,asymmetric = True):

	num  = np.asarray(num)
	sig  = np.sqrt(num + 1.0/4.) 	

	if asymmetric:
		return np.array([(sig +0.5)/num, (sig - 0.5)/num])
	else :
		return 1./np.sqrt(num)
    

def get_errorbars(fname  ):
    
    data  = np.loadtxt(fname)
    
    MininvhMsun = data[:,0]
    numcluster  = data[:,1]
    massfn      = data[:,2]
    calcB       = data[:,-2]
    calcM       = data[:,-1]
    
    return MininvhMsun , massfn, fracPoissonerrors(numcluster )*massfn, calcB, calcM

#MininvhMsun , massfn, massfnyerr = get_errorbars(fname)
#ax = plt.subplot()
def plotnormalized(fname , normals= None, normalnames=None, styles=None, usefile =  True):#, ps = None, cosmo = None):
    
    ax = plt.gca()
    #if axesobj != None:
    #       ax = axesobj
            
    MininvhMsun , massfn , massfnyerr,calcB, calcM = get_errorbars(fname)
    
    #print "YERR", massfnyerr
    #print np.shape(massfnyerr[0]), np.shape(massfnyerr[1]), np.shape( MininvhMsun), np.shape(massfn), np.shape(normals[0])
    #Mass = MininvhMsun*cosmo.h
    #normalB = hm.dndlnM(Mass , ps = ps , z= z, cosmo
    #print yerrba
    if normals == None:
        normval = 1.
        #ax.errorbar(MininvhMsun , massfn /normval , yerr = yerrbars/normval)
    if usefile:
        ax.errorbar( MininvhMsun, massfn/calcB , yerr = massfnyerr/calcB,  fmt  = "ks", label="Bhattacharya etal.")
        ax.errorbar( MininvhMsun, massfn/calcM , yerr = massfnyerr/calcM,  fmt  = "ro", label="MICE")
        
    else:
        for i, normval in enumerate(normals):
            name = normalnames [i] 
            stylename = styles[i]
            #print name, np.shape(normval), np.shape(massfn), np.shape(MininvhMsun), np.shape(normval), name
            ax.errorbar( MininvhMsun, massfn/normval , yerr = massfnyerr/normval, label = name, fmt  = stylename)
            #ax.errorbar( MininvhMsun, massfn/calcB , yerr = massfnyerr/calcB,  fmt  = "bd")


            #ax.errorbar( MininvhMsun, massfn/hm.dndlnM(M = masses, ps = psfromtkcbscb , z = 0.,bgtype="cb",cosmo= M000s), yerr = massfnyerr/normval, label = name, fmt  = stylename)
    ax.set_xscale('log')
    ax.grid(True)
    xvals = np.linspace(1.0e13,1.1e15,2)
    refval = 1.0
    bandvals = [-0.1,0.1]
    ax.plot(xvals, np.ones(len(xvals)),'k-',lw =2.0)
    ax.fill_between (xvals , refval + bandvals[0], refval + bandvals[1],color = 'gray', alpha =0.25)
    ax.set_xlim(1.0e13,1.0e15)
    ax.xaxis.set_ticklabels("",visible=False)
    ax.yaxis.set_ticks([0.9,1.0,1.1])
    #ax.yaxis.set_ticklabels("",visible=False)
    ax.set_ylim(0.85,1.15)
    ax.grid(True)

    return 0


# In[33]:

gs  = gridspec.GridSpec(3,1, height_ratios=[0.33,0.33,0.33], width_ratios=[1.,0.,0.])
ax0 = plt.subplot(gs[0])
fname = simdir+ "m000.499.fof_mf"
masses = get_errorbars(fname)[0]/M000s.h
#sigmaM_0 = psu.
dndlnm_B0 = hmf.dndlnM0(M = masses, ps = psfromtkcbscb , z = 0.,bgtype="cb",cosmo= M000s)/M000s.h**3
dndlnm_M0 = hmf.dndlnM0(M = masses, ps = psfromtkcbscb , z = 0.,bgtype="cb",cosmo = M000s, fittingform="MICE")/M000s.h**3 #*M000s.h**3

plotnormalized(fname, normals = [dndlnm_B0,dndlnm_M0],normalnames = ["Bhattacharya", "MICE"],styles =["ks","ro"])
#ax0.set_xscale('log')
#ax0.set_yscale('log')
ax1 = plt.subplot(gs[1])
psfromtkcbscb1 = psu.powerspectrum(koverh = None, pstype ="cb", sigma8type = "cb", asciifile =tkfile, cosmo = M000s,z = 1.0, method = "CAMBoutgrowth")

fname1 = simdir+ "m000.247.fof_mf"
masses1 = get_errorbars(fname1)[0]/M000s.h
dndlnm_B1 = hmf.dndlnM0(M = masses1, ps = psfromtkcbscb , z = 1.,bgtype="cb",cosmo= M000s,deltac= 1.684)/M000s.h**3/2.**3.
dndlnm_M1 = hmf.dndlnM0(M = masses1, ps = psfromtkcbscb , z = 1.,bgtype="cb",cosmo = M000s, fittingform="MICE")/M000s.h**3./2**3.0 #*M000s.h**3

#dndlnm_B1 = hmf.dndlnM(M = masses1, ps = psfromtkcbscb1 , z = 1.,bgtype="cb",cosmo= M000s,deltac= 1.684)/M000s.h**3
#dndlnm_M1 = hmf.dndlnM(M = masses1, ps = psfromtkcbscb1 , z = 1.,bgtype="cb",cosmo = M000s, fittingform="MICE")/M000s.h**3. #*M000s.h**3
plotnormalized(fname1, normals = [dndlnm_B1,dndlnm_M1],normalnames = ["Bhattacharya", "MICE"],styles =["ks","ro"])

#plotnormalized(simdir+ "m000.247.fof_mf")
ax2 = plt.subplot(gs[2])
psfromtkcbscb2 = psu.powerspectrum(koverh = None, pstype ="cb", sigma8type = "cb", asciifile =tkfile, cosmo = M000s,z = 2.0, method = "CAMBoutgrowth")
fname2 = simdir+ "m000.163.fof_mf"
masses2 = get_errorbars(fname2)[0]/M000s.h

#print "Masses", masses2
#print "Eavltest", np.shape(masses2)
dndlnm_B2 = hmf.dndlnM0(M = masses2, ps = psfromtkcbscb , z = 2.,bgtype="cb",cosmo= M000s, deltac=1.686)/M000s.h**3/3**3.0
dndlnm_M2 = hmf.dndlnM0(M = masses2, ps = psfromtkcbscb , z = 2.,bgtype="cb",cosmo = M000s, fittingform="MICE")/M000s.h**3/3.0**3 #*M000s.h**3

#dndlnm_B2 = hmf.dndlnM(M = masses2, ps = psfromtkcbscb2 , z = 2.,bgtype="cb",cosmo= M000s, deltac=1.686)/M000s.h**3/27.0
#dndlnm_M2 = hmf.dndlnM(M = masses2, ps = psfromtkcbscb2 , z = 2.,bgtype="cb",cosmo = M000s, fittingform="MICE")/M000s.h**3/27. #*M000s.h**3
#print "Hello what is the problem?"
#print np.shape(masses2), np.shape(psfromtkcbscb2), np.shape(dndlnm_B2),np.shape(dndlnm_M2)
#plotnormalized(fname2, normals = [dndlnm_B1,dndlnm_M1],normalnames = ["Bhattacharya", "MICE"],styles =["ks","ro"])

#print "Bye"
#print "Eavl"#, np.shape(masses2) #,np.shape(psfromtkcbscb2), np.shape(dndlnm_B2) , np.shape(dndlnm_M2)

#plotnormalized(fname2, normals = [dndlnm_B2,dndlnm_M2],normalnames = ["Bhattacharya", "MICE"],styles =["ks","ro"])

plotnormalized(fname2, normals = [dndlnm_B2,dndlnm_M2],normalnames = ["Bhattacharya", "MICE"],styles =["ks","ro"])

#plotnormalized(simdir + "m000.163.fof_mf")
ax2.set_xscale('log')
ax2.set_xlabel(r'Mass ($h^{-1}M_\odot$)' )
ax1.set_ylabel("ratio of massfn")
ax2.legend(loc='upper right', numpoints =1)
ax1.text(1.0, 0.0,"z = 2.")
#axt2 = ax2.twiny()
#axt2.set_xlabel(r"z=2")
#axt1 = ax1.twiny()
#axt1.set_xlabel(r"z=1")
#axt0 = ax0.twiny()
#axt0.set_xlabel(r"z=0")

plt.tight_layout(h_pad = -0.1, w_pad = -0.1)
#plt.savefig("M000_Formula.pdf")
plt.savefig("M000_SumanFormula.pdf")


# Out[33]:

#     0.333 0.788 0.807 1.795
#     0.308553824609 0.782556886404 0.807 1.795
#     0.295094350742 0.779390315288 0.807 1.795
# 

# image file:

# In[17]:

plt.plot(masses,1.0/hmf.dndlnM0(M = masses, ps = psfromtkcbscb , z = 2.,bgtype="cb",cosmo= M000s, deltac=1.686)*hmf.dndlnM(M = masses, ps = psfromtkcbscb2, z= 2.0, bgtype = "cb", cosmo= M000s, deltac=1.686),"ks")
plt.xscale('log')


# Out[17]:

#     0.295094350742 0.779390315288 0.807 1.795
#     0.295094350742 0.779390315288 0.807 1.795
# 

# image file:

# In[ ]:

dndlnm_B1c = dndlnm_B1*0.12
plotnormalized(fname1, normals = [dndlnm_B1c,dndlnm_M1],normalnames = ["Bhattacharya", "MICE"],styles =["ks","ro"])


print dndlnm_B1 , dndlnm_B1c


# In[38]:


dndlnm_B2 = hm.dndlnM(M = masses2, ps = psfromtkcbscb2 , z = 2.,bgtype="cb",cosmo= M000s, deltac=1.686)/M000s.h**3
dndlnm_M2 = hm.dndlnM(M = masses2, ps = psfromtkcbscb2 , z = 2.,bgtype="cb",cosmo = M000s, fittingform="MICE")/M000s.h**3
print dndlnm_B2, dndlnm_M2
print np.shape(dndlnm_B2)
print 1./8.


# Out[38]:

#     0.295094350742 0.779390315288 0.807 1.795
#     [  8.69840997e-03   5.74194925e-03   3.62316898e-03   2.18669780e-03
#        1.27152405e-03   7.13534935e-04   3.82356121e-04   1.92706374e-04
#        9.00314410e-05   3.81860774e-05   1.49713792e-05   5.08301266e-06
#        1.40307142e-06   3.45923170e-07   6.94611686e-08   1.67266650e-08] [  8.09230756e-03   5.33340159e-03   3.36490446e-03   2.03425732e-03
#        1.18760520e-03   6.71018535e-04   3.63331163e-04   1.85832054e-04
#        8.85716560e-05   3.85895834e-05   1.56552636e-05   5.55795206e-06
#        1.62826867e-06   4.31161409e-07   9.46077053e-08   2.47736866e-08]
#     (16,)
#     0.125
# 

# In[109]:

plt.plot(masses, hm.dndlnM(M = masses, ps = psfromtkcbscb, cosmo = M000s)/psu.dndlnM(masses, ps = psfromtkcbscb,cosmo = M000))
plt.xscale('log')
plt.ylim(0.9,1.1)
print M000s.h**3

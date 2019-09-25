from minirats.utils.py.cellutils import py_cell_utils as cu
from astro import readRamses as rr
from astro import constants as c
from astro import SEDutils as su
from astro import cloudy as cl
from algebra import utils as alg
import numpy as np
import math
import os
import random


def compute_intensity_Cloudy(sedDir, ramses_intensity):

	ryd, Fnu = np.genfromtxt(sedDir, unpack=True)	
	nu = c.ryd_to_Hz*ryd

	return 1/np.trapz(Fnu/(c.h_cgs*nu), nu) * np.trapz(Fnu, nu)
	


ramDir = '/home/evol/mauerhof/RamsesFiles/idealized/G8/newRamses/Dust/'
timestep = 12
center = [0.5, 0.5, 0.5]
radius = 0.06
sedDir = './seds/'

indices = [1, 5, 6, 7, 8, 9] #Density, pressure, metallicity, xHII, xHeII, xHeIII


xH = rr.read_xH(ramDir, timestep)
yHe = 1. - xH
units = rr.read_units(ramDir, timestep)
groups = rr.read_groups(ramDir, timestep)
lmax = rr.read_lmax(ramDir, timestep)
nBins = len(groups)
binIndices = range(10,10+4*nBins,4)


ncells = cu.count_cells(ramDir,timestep,lmax,center,radius)


cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,indices+binIndices,center,radius,True)

cells0,cell_pos0,cell_l0 = cu.read_cells_hydro(ramDir+'restart4/testLoadBalancing/',timestep+1,lmax,ncells,[10],center,radius,True)


ramses_to_cloudy = np.array([])
for i in range(4):
	ryd, Fnu = np.genfromtxt(sedDir+'sed_cloudy'+str(i)+'.sed', unpack=True)
	nu = c.ryd_to_Hz*ryd
	if (i > 0):
		ramses_to_cloudy = np.append(ramses_to_cloudy, 1/np.trapz(Fnu/(c.h_cgs*nu), nu) * np.trapz(Fnu, nu))
	else:
		ramses_to_cloudy = np.append(ramses_to_cloudy, 1/np.trapz(Fnu[nu > c.ryd_to_Hz*0.599]/(c.h_cgs*nu[nu > c.ryd_to_Hz*0.599]), nu[nu > c.ryd_to_Hz*0.599]) * np.trapz(Fnu[nu > c.ryd_to_Hz*0.599], nu[nu > c.ryd_to_Hz*0.599]))
	
nH = xH*cells[:,0]#*units[1]/c.mp_cgs
xHII = cells[:,3]
xHeII = cells[:,4]
xHeIII = cells[:,5]
			
Tgas = cells[:,1]/cells[:,0]*1.67e-24/1.38062e-16*(units[0]/units[2])**2 / (xH*(1+xHII[:]) + 0.25*yHe*(1+xHeII[:]+2*xHeIII[:]))     #mp/kb * (unit_length/unit_time)^2  is the conversion factor
#if (Tgas < 1d6):
intensities = np.zeros((4, ncells))
intensities[0,:] = ramses_to_cloudy[0]*cells0[:,0]*units[0]/units[2]

for i in range(1,4):
	intensities[i,:] = ramses_to_cloudy[i]*cells[:,len(indices)+i-1]*units[0]/units[2]

np.savetxt('input.txt', np.transpose([nH, Tgas, cells[:,2], intensities[0,:], intensities[1,:], intensities[2,:], intensities[3,:]]), ['%.5e', '%.5e', '%.5e', '%.5e', '%.5e', '%.5e', '%.5e'])


	







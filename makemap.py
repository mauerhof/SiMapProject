from minirats.utils.py.cellutils import py_cell_utils as cu
from astro import readRamses as rr
from matplotlib import pyplot as plt
import numpy as np
import math


# read cellmax = 0 # read all
lmax = 13
indices = [1,5,6,7,8,9] #Density, pressure, metallicity, xHII, xHeII, xHeIII
 
center = [0.5,0.5,0.5]
radius = 1                     

#ramDir = '/Users/mauerhof/Documents/RamsesFiles/idealized/G8/lmax13Ngroups7/'
ramDir = './'
xlim = [0.47,0.53]
ylim = [0.47,0.53]
zlim = [0.495,0.505]

sigmaSi = np.genfromtxt('SiCSN.dat', unpack=True)



for timestep in range(13,16):  


	units = rr.read_units(ramDir, timestep)   #length, d (?), time
	#print units

	groups = rr.read_groups(ramDir, timestep)
	#print groups


	nBins = len(groups)
	binIndicies = range(10,10+4*nBins,4)
		      

	ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
	cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,indices+binIndicies,center,radius,True)

	Tgas = cells[:,1]/cells[:,0]*1.67e-24/1.38062e-16*(units[0]/units[2])**2
	
	rateSi = np.zeros([ncells,4])
	for i in range(ncells):
		for j in range(4):
			rateSi[i,j] = np.sum(cells[i,6:6+nBins:1]*sigmaSi[j,:])
	rateSi = rateSi*units[0]/units[2]
	romanNum = ['I', 'II', 'III', 'IV']
	
	normalizedMass = cells[:,0]*np.power(8*np.ones_like(cells[:,0]),-cell_l)
	normalizedMass = normalizedMass/np.sum(normalizedMass)


	nx, ny = cu.get_map_nxny(lmax, xlim[0], xlim[1], ylim[0], ylim[1])
	
	for i in range(4):
		TempMap, w = cu.make_map(lmax, False, xlim[0], xlim[1], ylim[0], ylim[1], zlim[0], zlim[1], np.log10(rateSi[:,i]), cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
		plt.imshow(TempMap, origin='lower')
		#plt.clim(-12.2, -6)
		cbar = plt.colorbar()
		cbar.set_label(r'Photoionization rate of Si'+romanNum[i]+', log[$s^{-1}$]')
		plt.savefig('Si'+romanNum[i]+'rate_Dust_'+str(timestep)+'.png')
		plt.close()

	TempMap, w = cu.make_map(lmax, False, xlim[0], xlim[1], ylim[0], ylim[1], zlim[0], zlim[1], np.log10(Tgas), cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
	plt.imshow(TempMap, origin='lower')
	plt.colorbar()
	plt.savefig('T_Map_'+str(timestep)+'.png') 
	plt.close()

	TempMap, w = cu.make_map(lmax, False, xlim[0], xlim[1], ylim[0], ylim[1], zlim[0], zlim[1], cells[:,3], cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
	plt.imshow(TempMap, origin='lower')
	plt.colorbar()
	plt.savefig('HII_Map_'+str(timestep)+'.png') 
	plt.close()

	TempMap, w = cu.make_map(lmax, False, xlim[0], xlim[1], ylim[0], ylim[1], zlim[0], zlim[1], cells[:,2], cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
	plt.imshow(TempMap, origin='lower')
	plt.colorbar()
	plt.savefig('Z_Map_'+str(timestep)+'.png') 
	plt.close()
	
	plt.figure(figsize=(14, 8))
	plt.scatter(np.log10(cells[:,0]), np.log10(Tgas), s=0.05, marker='.', c=np.log10(normalizedMass))
	cbar = plt.colorbar()
	cbar.set_label('f (mass)')
	plt.xlabel('log(nH[cm-3])')
	plt.ylabel('log(T [K])')
	plt.savefig('phaseDiagram_'+str(timestep)+'.png')
	plt.close()

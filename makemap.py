from minirats.utils.py.cellutils import py_cell_utils as cu
from astro import readRamses as rr
from astro import constants as c
import matplotlib
# ~ matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import math

#Directory where Ramses output is
ramDir = '/Users/mauerhof/Documents/RamsesFiles/idealized/G8/newRamses/Dust/'

#Please set the indices here
indices = [1,5,6,7,8,9] #Density, pressure, metallicity, xHII, xHeII, xHeIII
 
#Center and radius of the region extracted by minirats py_cell_utils 
center = [0.5,0.5,0.5]
radius = 6e-2                   

#Region for the plot
xlim = [0.47,0.53]
ylim = [0.47,0.53]
zlim = [0.4999,0.5001]  #Thin slice




for timestep in range(12,13):  


	units = rr.read_units(ramDir, timestep)   #length, density, time
	#print units

	groups = rr.read_groups(ramDir, timestep)
	#print groups
	
	lmax = rr.read_lmax(ramDir, timestep)
	boxlen = rr.read_boxlen(ramDir, timestep)
	xH = rr.read_xH(ramDir, timestep)


	nBins = len(groups)
	binIndicies = range(max(indices)+1,max(indices)+1+4*nBins,4)
	

	ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
	print (ncells)
	cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,indices+binIndicies,center,radius,True)

	Tgas = cells[:,1]/cells[:,0]*1.67e-24/1.38062e-16*(units[0]/units[2])**2 
	
	nH = xH*cells[:,0]*units[1]/c.mp_cgs


	nx, ny = cu.get_map_nxny(lmax, xlim[0], xlim[1], ylim[0], ylim[1])

	TempMap, w = cu.make_map(lmax, False, xlim[0], xlim[1], ylim[0], ylim[1], zlim[0], zlim[1], np.log10(nH), cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
	plt.style.use('dark_background')
	plt.imshow(TempMap, origin='lower', extent=[-(xlim[1]-xlim[0])*boxlen/2, (xlim[1]-xlim[0])*boxlen/2, -(ylim[1]-ylim[0])*boxlen/2, (ylim[1]-ylim[0])*boxlen/2], cmap=plt.get_cmap('seismic'))
	# ~ plt.clim(-20, -5)
	plt.xlabel('[kpc]')
	plt.ylabel('[kpc]')
	cbar = plt.colorbar()
	cbar.set_label(r'$\mathrm{log}(n_H [\mathrm{cm}^{-3}])$', fontsize=20)
	# ~ plt.savefig('./Si'+romanNum[i]+'_density.png')
	plt.close()

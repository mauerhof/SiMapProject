from minirats.utils.py.cellutils import py_cell_utils as cu
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import math
from astro import constants as c
from astro import utils as u

# ~ pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
# ~ matplotlib.rcParams.update(pgf_with_rc_fonts)

bins = [8.150, 13.600, 16.350, 24.590, 33.490, 45.140, 54.420, 1e10]
nBins = len(bins)-1
unit_length = 0.308567758128200e22
unit_time = 0.470430312423675e15

lmax = 13
hydroIndices = [1,5,6,7,8,9] #density, pressure, metallicity, ixHII, ixHeII, ixHeIII
binIndicies = range(10,10+4*nBins,4)
indices = hydroIndices + binIndicies

 
center = [0.5,0.5,0.5]
radius = 1                     

ramDir = '.'

timestep = 12                   

ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,indices,center,radius,True)

normalizedMass = cells[:,0]*np.power(8*np.ones_like(cells[:,0]),-cell_l)
normalizedMass = normalizedMass/np.sum(normalizedMass)


Galaxycells = cells[(cell_pos[:,0]-0.5)**2 + (cell_pos[:,1]-0.5)**2 + (cell_pos[:,2]-0.5)**2 < (255./8192.)**2]
Galaxycell_pos = cell_pos[(cell_pos[:,0]-0.5)**2 + (cell_pos[:,1]-0.5)**2 + (cell_pos[:,2]-0.5)**2 < (255./8192.)**2]
Galaxycell_l = cell_l[(cell_pos[:,0]-0.5)**2 + (cell_pos[:,1]-0.5)**2 + (cell_pos[:,2]-0.5)**2 < (255./8192.)**2]

GalaxynormalizedMass = Galaxycells[:,0]*np.power(8*np.ones_like(Galaxycells[:,0]),-Galaxycell_l)
GalaxynormalizedMass = GalaxynormalizedMass/np.sum(GalaxynormalizedMass)

# ~ print ncells
# ~ print Galaxycells[:,0].size
# ~ print ncells
# ~ print cells[:,0].min()    
# ~ print cells[:,0].max()
# ~ print cells[:,1].min()
# ~ print cells[:,1].max()
# ~ print cells[:,2].min()
# ~ print cells[:,2].max()
# ~ print cells[:,3].min()
# ~ print cells[:,3].max()
# ~ print cells[:,4].min()
# ~ print cells[:,4].max()
# ~ print cell_pos[1:100:1,0] 
# ~ print cell_pos[:,0].max()
# ~ print cell_l.max() 


Tgas = cells[:,1]/cells[:,0]*1.67e-24/1.38062e-16*(unit_length/unit_time)**2
GalaxyTgas = Galaxycells[:,1]/Galaxycells[:,0]*1.67e-24/1.38062e-16*(unit_length/unit_time)**2

ssp = u.readRamsesSEDs('/Users/mauerhof/Documents/seds/bpass100')
flux = ssp['spectra'][:,9,1]
wavelength = ssp['lambdaBins']

sigmaSiII = np.empty(nBins)
for i in range(nBins):
	sigmaSiII[i] = u.csn(flux, wavelength, c.hc_AeV/bins[i], c.hc_AeV/bins[i+1], 'Si', 2)

rate = np.zeros(ncells)
for i in range(ncells):
	rate[i] = np.sum(cells[i,6:6+nBins:1]*sigmaSiII[:])
rate = rate*unit_length/unit_time
print cells[1:100:1, 6:13:1]*unit_length/unit_time


# ~ plt.rc('text', usetex=True)
# ~ plt.rc('font', family='serif')

xylim = [0.46,0.54]

nx, ny = cu.get_map_nxny(lmax, xylim[0], xylim[1], xylim[0], xylim[1])
weight = np.ones(ncells)
TempMap, w = cu.make_map(lmax, False, xylim[0], xylim[1], xylim[0], xylim[1], 0.495, 0.505, np.log10(rate), weight, cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
plt.imshow(TempMap, origin='lower')
cbar = plt.colorbar()
cbar.set_label(r'Photoionization rate of SiII log[$s^{-1}$]')
#plt.savefig('photoRateSiII_Ngroups7.png') 
#plt.show()




# plt.figure(figsize=(14, 8))
# plt.scatter(np.log10(Galaxycells[:,0]), np.log10(GalaxyTgas), s=0.05, marker='.', c=np.log10(GalaxynormalizedMass))
# cbar = plt.colorbar()
# cbar.set_label('f (mass)')
# plt.xlabel('log(nH[cm-3])')
# plt.ylabel('log(T/mu [K])')
# plt.savefig('phaseDiagramGalaxy.png')
#plt.show()
 

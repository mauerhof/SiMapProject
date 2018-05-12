from minirats.utils.py.cellutils import py_cell_utils as cu
from astro import readRamses as rr
from astro import constants as c
from astro import SEDutils as su
from algebra import utils as alg
import numpy as np
import math
from mpi4py import MPI


comm = MPI.COMM_WORLD
size = comm.Get_size()
print size
rank = comm.Get_rank()


#  This part contains variables that the user chooses for each output.  Another thing to choose at cloudy S definition
# ~ ramDir = '/Users/mauerhof/Documents/RamsesFiles/idealized/G8/lmax13Ngroups3/'
ramDir = './'
timestep = 16
center = [0.5, 0.5, 0.5]
radius = 1
ions = 4

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

# ~ cloudyE = groups[:-1] + np.diff(groups)/2
# ~ cloudyE = np.append(cloudyE, groups[-1]+5)
cloudyE = np.array([9., 17., 34., 46.])

cloudyC = np.empty([ions,ions])   #First index ion species (SiI, SII, ...),  second index for bin
for i in range(ions):
	for j in range(ions):
		cloudyC[i,j] = su.photoCS(c.hc_AeV/cloudyE[j], 'Si', i+1)/cloudyE[j]/c.ev_to_erg
		
ramsesCSN = np.genfromtxt('SiCSN.dat', unpack=True)

ratesRamses = np.empty([ions])
for icells in range(1):
	for i in range(ions):
		ratesRamses[i] = np.sum(cells[icells,len(indices):len(indices)+nBins:1]*ramsesCSN[i,:])*units[0]/units[2]
	
	cloudyI = alg.solveTriangularMatrix(ratesRamses, cloudyC)
	if rank==0: 
		print rank, cloudyI

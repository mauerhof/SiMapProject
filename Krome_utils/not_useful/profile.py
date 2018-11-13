from minirats.utils.py.cellutils import py_cell_utils as cu
from matplotlib import pyplot as plt
import numpy as np
import math
from astro import constants as c


hc_AeV = 12398.4  #A * eV
hc_cgs = 1.98644e-16
c_cgs = 2.99792e10
h_cgs = hc_cgs/c_cgs
HI, HeI, HeII = 0, 1, 2
Eion = np.array([13.6, 24.59, 54.42])
lambdaIon = np.array([911.75, 504.257, 227.837])
#nuIon = [c_cgs/lambdaIon_A[HI]/1e-8, c_cgs/lambdaIon_A[HeI]/1e-8, c_cgs/lambdaIon_A[HeII]/1e-8]
erg_to_ev = 6.2415e11
Lsun_cgs = 3.84e33
pc_to_cm = 3.08568e18
kb_cgs = 1.38065e-16



def photoCrossSection(lambdaBins):    #From Verner et al. (1996),  H and He
	
	epsilon0 = [0.4298, 13.61, 1.720] #eV
	sigma0 = [5.475e-14, 9.492e-16, 1.369e-14] #cm^2
	P = [2.963, 3.188, 2.963]
	ya = [32.88, 1.469, 32.88]
	yw = [0, 2.039, 0]
	y0 = [0, 0.4434, 0]
	y1 = [0, 2.136, 0]
	
	
	n = len(lambdaBins)
	
	xx = np.empty([3, n])
	yy = np.empty([3, n])
	sigma = np.empty([3, n])
	for element in range(3):
		for i in range(n):
			if lambdaBins[i] > lambdaIon[element]:
				xx[element, i] = 0
				yy[element, i] = 0
				sigma[element, i] = 0
			else:
				xx[element, i] = hc_AeV/lambdaBins[i]/epsilon0[element] - y0[element]
				yy[element, i] = math.sqrt(xx[element, i]**2 + y1[element]**2)
				sigma[element, i] = sigma0[element]*((xx[element, i] - 1)**2 + yw[element]**2) * math.pow(yy[element, i], 0.5*P[element]-5.5)/math.pow(1 + math.sqrt(yy[element,i]/ya[element]), P[element])
		
	
	return sigma




def Si():
	SiIon = c.SiIon[c.SiIon < np.amax(c.SiIon)]
	partition = c.HIon
	partition = np.append(partition, c.HeIon)
	partition = np.append(partition, SiIon)
	partition = c.hc_AeV/partition
	partition = np.append(partition, np.array([0]))
	#partition = np.append(partition, np.array([1215]))
	#partition = np.append(partition, np.array([630]))
	partition = partition[np.argsort(partition)]
	
	return partition




if __name__ == '__main__':


	# read cells
	lmax = 6 # read all levels 

	idens,ipre,ixhii,ixheii,ixheiii = 1,6,7,8,9     #Density, pressure, fraction of HII, HeII, HeIII
	center = [0.5,0.5,0.5]
	radius = 1.0

	# RT fields
	iNp1,iNp2,iNp3,iNp4,iNp5,iNp6,iNp7,iNp8,iNp9 = 10,14,18,22,26,30,34,38,42
	readRT=True

	ramDir = '/Users/mauerhof/Documents/RamsesFiles/stromgren/'
	timestep = 17  #snapshot   
	ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
	cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,[idens,ipre,ixhii,ixheii,ixheiii,iNp1,iNp2,iNp3,iNp4,iNp5,iNp6,iNp7,iNp8,iNp9],center,radius,readRT)
	#print (ncells)
	
	
	rank,E,HI,HeI,SiI,HII,HeII,SiII,HeIII,SiIII,SiIV,SiV,SiVI,SiVII = np.genfromtxt('./outstromgrenSi.dat', unpack=True)
	
	

	ymin = cell_pos[:,1].min()
	zmin = cell_pos[:,2].min()

	X = cell_pos[:,0][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]    #Extract the list of x axis points,  for y and z at 0
	# #print X.max()
	radius = X*3.086e21 #converts to cm
	# print cells[:,5]
	SiI_X = SiI[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	SiII_X = SiII[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	SiIII_X = SiIII[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	SiIV_X = SiIV[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	SiV_X = SiV[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	SiVI_X = SiVI[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	SiVII_X = SiVII[(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]


	data = np.column_stack((radius, SiI_X, SiII_X, SiIII_X, SiIV_X, SiV_X, SiVI_X, SiVII_X))
	data = data[np.argsort(data[:, 0])]  #Rearrange the data to have increasing radius
	
	Sifactor = 1./3.467e-6
	
	fig1, ax1 = plt.subplots(figsize=(9,7))
	ax1.plot(data[:,0], Sifactor*data[:,1], 'r', label='SiI Krome')
	ax1.plot(data[:,0], Sifactor*data[:,2], 'b', label='SiII Krome')
	ax1.plot(data[:,0], Sifactor*data[:,3], 'g', label='SiIII Krome')
	ax1.plot(data[:,0], Sifactor*data[:,4], 'orange', label='SiIV Krome')
	ax1.plot(data[:,0], Sifactor*data[:,5], 'brown', label='SiV Krome')
	ax1.plot(data[:,0], Sifactor*data[:,6], 'purple', label='SiVI Krome')
	ax1.set(xlabel='Radius [cm]', ylabel='Silicone fractions', title='Stromgren sphere, Bpass100_z0.002_6.3Myr, 1e50 photons/s, 1e4 K')

	legend=ax1.legend(loc='center right', prop={'size': 16})
	ax1.grid()

	plt.show()
	#fig1.savefig('./Si_7bins_Krome.png')
	
	


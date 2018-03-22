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


def readRamsesSEDs(sedDir):
   """Read SED in ramses format and return

   Parameters:
   ----------------------------------------------------------------------
   sedDir: Directory containing the SED tables
   """
   from scipy.io import FortranFile
   # Read metallicity bins
   ZFile = open(sedDir+'/metallicity_bins.dat', 'r')
   nZ = eval(ZFile.readline())
   ZBins = []
   for Z in range(0,nZ): 
	   ZBins.append(eval(ZFile.readline()))
   ZFile.close()

   # Read age bins
   ageFile = open(sedDir+'/age_bins.dat', 'r')
   nAge = eval(ageFile.readline())
   ageBins = []
   for age in range(0,nAge): ageBins.append(eval(ageFile.readline()))
   ageFile.close()

   # Read wavelength bins and spectra
   sedFile = FortranFile(sedDir+'/all_seds.dat','r')
   nLambda = sedFile.read_ints()[0]
   lambdaBins = sedFile.read_reals()
   spectra = np.empty([nLambda,nAge,nZ]) 
   for iZ in range(0,nZ):
       for iAge in range(0,nAge):
           spectrum = sedFile.read_reals()
           spectra[:,iAge,iZ] = Lsun_cgs*spectrum  

   return {'ZBins':ZBins, 'ageBins':ageBins, 'lambdaBins':lambdaBins,'spectra':spectra, 'nLambda':nLambda}   #spectra are in erg/s/A/Msun



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



def csn(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 250, 500, 911.75) -> 3 bins
	
	sigma = photoCrossSection(wavelength)
	csn_array1 = flux*sigma*wavelength*1e-8/hc_cgs
	csn_array2 = flux*wavelength*1e-8/hc_cgs
	
	csn = np.empty([len(partition)-1, 3])
	for bins in range(len(partition)-1):
		for element in range(3):
			#print wavelength[(partition[bins] < wavelength) & (wavelength < partition[bins+1])]
			csn[-bins-1][element] = np.trapz(csn_array1[element][(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])/(np.trapz(csn_array2[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])]))


	return csn
	
	
def cse(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 250, 500, 911.75) -> 3 bins
	
	sigma = photoCrossSection(wavelength)
	cse_array1 = flux*sigma
	cse_array2 = flux
	
	cse = np.empty([len(partition)-1, 3])
	for bins in range(len(partition)-1):
		for element in range(3):
			#print wavelength[(partition[bins] < wavelength) & (wavelength < partition[bins+1])]
			cse[-bins-1][element] = np.trapz(cse_array1[element][(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])/(
			np.trapz(cse_array2[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])]))


	return cse



def egy(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 227.84, 504.26, 911.75) -> 3 bins
	
	egy_array1 = flux
	egy_array2 = flux*wavelength*1e-8/hc_cgs
	
	egy = np.empty([len(partition)-1])
	for bins in range(len(partition)-1):
		#print wavelength[(partition[bins] < wavelength) & (wavelength < partition[bins+1])]
		egy[bins] = np.trapz(egy_array1[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])/(
		np.trapz(egy_array2[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])]))


	return np.flip(egy,0)*erg_to_ev #in eV
	
	
def photons(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 227.84, 504.26, 911.75) -> 3 bins
	
	#sigma = photoCrossSection(wavelength)
	photons_array = flux*wavelength*1e-8/hc_cgs
	
	photons = np.empty([len(partition)-1])
	for bins in range(len(partition)-1):
		#print wavelength[(partition[bins] <= wavelength) & (wavelength < partition[bins+1])]
		photons[bins] = np.trapz(photons_array[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])


	return photons #in eV



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

	ramDir = '/Users/mauerhof/Documents/Ramses/stromgren/'
	timestep = 17  #snapshot   
	ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
	cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,[idens,ipre,ixhii,ixheii,ixheiii,iNp1,iNp2,iNp3,iNp4,iNp5,iNp6,iNp7,iNp8,iNp9],center,radius,readRT)
	#print (ncells)
	

	ymin = cell_pos[:,1].min()
	zmin = cell_pos[:,2].min()
	#print cells[:,0]
	#print cells[:,1]
	#print cells[:,1]/cells[:,0]*1.66e-24/1.38062e-16*(3.09e21/3.1556926e13)**2
	#print 1.66e-24/1.38062e-16*(3.09e21/3.1556926e13)**2
	X = cell_pos[:,0][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]    #Extract the list of x axis points,  for y and z at 0
	#print X.max()
	radius = X*3.086e21 #converts to cm
	print cells[:,5]
	
	HII = cells[:,2][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]  #Those lines take the ionization fractions corresponding to the x axis points X
	HI = 1-HII

	HeII = cells[:,3][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	HeIII = cells[:,4][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	HeI = 1-HeII-HeIII

	N1 = cells[:,5][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]*3.09e21/c.myr_to_s  #c*N for each photon Bin,  still y and z at 0
	N2 = cells[:,6][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N3 = cells[:,7][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N4 = cells[:,8][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N5 = cells[:,9][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N6 = cells[:,10][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N7 = cells[:,11][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N8 = cells[:,12][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	N9 = cells[:,13][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]
	

	data = np.column_stack((radius, HI, HII, HeI, HeII, HeIII, N1, N2, N3, N4, N5, N6, N7, N8, N9))
	data = data[np.argsort(data[:, 0])]  #Rearrange the data to have increasing radius
	
	
	#*******************************
	#To compute the optical depth
	Hdensity = 0.1
	Hedensity = 0.0789*Hdensity
	
	
	ssp = readRamsesSEDs('/Users/mauerhof/Documents/seds/bpass100')
	wavelength = ssp['lambdaBins']
	
	partition = Si()
	
	photonsRamses = 1e50
	flux = ssp['spectra'][:, 9, 1]
	
	Bins = partition.size - 1
	

	sigma = csn(flux, wavelength, partition)
	energies = egy(flux, wavelength, partition)

	tau = np.empty([Bins, data[:,0].size])
	for i in range(data[:,0].size):
		for j in range(Bins):
			tau[j][i] = sigma[j][0]*np.trapz(Hdensity*data[:,1][0:i+1:1], data[:,0][0:i+1:1]) + sigma[j][1]*np.trapz(Hedensity*data[:,3][0:i+1:1], data[:,0][0:i+1:1]) + sigma[j][2]*np.trapz(Hedensity*data[:,4][0:i+1:1], data[:,0][0:i+1:1])

	probAbsorption = np.zeros_like(tau)
	for i in range(Bins):
		probAbsorption[i,:] = 1 - np.exp(-tau[i,:])
	#probAbsorption = 1 - np.exp(-tau)

	probAbsorption_mesured = np.zeros_like(probAbsorption)
	for i in range(Bins):
		probAbsorption_mesured[i,:] = 1 - data[:,6+i]*data[:,0]**2/np.max(data[:,6+i]*data[:,0]**2)

	#************************************





	# fig2, ax2 = plt.subplots(figsize=(9,7))
	# #ax2.plot(data[:,0], data[:,9]*data[:,0]**2, 'b', label='c_r * Nphotons*r^2')
	# #ax2.plot(data[:,0], 1-data[:,1], 'b', label='HI')
	# ax2.plot(data[:,0], data[:,2], 'b--', label='HII')
	# #ax2.plot(data[:,0], 1-data[:,2]-data[:,3], 'r', label='HeI')
	# ax2.plot(data[:,0], data[:,4], 'r--', label='HeII')
	# ax2.plot(data[:,0], data[:,5], 'r:', label='HeIII')
	# ax2.set(xlabel='Radius [cm]', ylabel='Ionization fractions', title='Stromgren sphere, 3Bins, 1e49 ionizing photons/s')
	# # ax.set_ylim(1.0e-10)
	# # #ax.set_xlim(3.0, 4.0)
	# legend=ax2.legend(loc='center right', prop={'size': 16})
	# ax2.grid()
	# #fig2.savefig('ionization_fractions.png')


	# fig4, ax4 = plt.subplots(figsize=(9,7))
	# for i in range(9):
		# ax4.plot(data[:,0], probAbsorption[i]-probAbsorption_mesured[i], label=str(round(energies[i],1)) +' eV')
		# #ax4.plot(data[:,0], probAbsorption_mesured[i], 'r', label='1 - n(r)/nmax, 18 eV')
	# # ax3.plot(data[:,0], data[:,8], 'r', label='22 eV')
	# # ax3.plot(data[:,0], data[:,12], 'g', label='60 eV')
	# ax4.set(xlabel='Radius [cm]', ylabel='probability of absorption', title='probability of absorption')
	# # ax.set_ylim(1.0e-10)
	# # #ax.set_xlim(3.0, 4.0)
	# legend=ax4.legend(loc='lower right', prop={'size': 14})
	# ax4.grid()
	# #fig3.savefig('photonsPerCC.png')


	# fig3, ax3 = plt.subplots(figsize=(9,7))
	# ax3.plot(data[:,0], data[:,4], 'b', label='9 eV')
	# ax3.plot(data[:,0], data[:,8], 'r', label='22 eV')
	# ax3.plot(data[:,0], data[:,12], 'g', label='60 eV')
	# ax3.set(xlabel='Radius [cm]', ylabel='c * Nphotons', title='Photon propagation')
	# # ax.set_ylim(1.0e-10)
	# # #ax.set_xlim(3.0, 4.0)
	# legend=ax3.legend(loc='center right', prop={'size': 16})
	# ax3.grid()
	# #fig3.savefig('photonsPerCC.png')


	
	fig1, ax1 = plt.subplots(figsize=(9,7))
	#for i in range(4,13,1):
	#ax1.plot(data[:,0], data[:,i]*data[:,0]**2)
	ax1.plot(data[:,0], data[:,6]*np.power(data[:,0],2)*4*math.pi, 'b', label='9 eV')  #data[:,6] is c*N for the 9eV photons
	# ax1.plot(data[:,0], data[:,5]/data[0,5]*data[:,0]**2/data[0,0]**2, 'r', label='11 eV')
	# ax1.plot(data[:,0], data[:,6]/data[0,6]*data[:,0]**2/data[0,0]**2, 'g', label='15 eV')
	# ax1.plot(data[:,0], data[:,7]/data[0,7]*data[:,0]**2/data[0,0]**2, 'violet', label='18 eV')
	# ax1.plot(data[:,0], data[:,8]/data[0,8]*data[:,0]**2/data[0,0]**2, 'orange', label='22 eV')
	# ax1.plot(data[:,0], data[:,9]/data[0,9]*data[:,0]**2/data[0,0]**2, 'pink', label='28 eV')
	# ax1.plot(data[:,0], data[:,10]/data[0,10]*data[:,0]**2/data[0,0]**2, 'maroon', label='38 eV')
	# ax1.plot(data[:,0], data[:,11]/data[0,11]*data[:,0]**2/data[0,0]**2, 'navy', label='48 eV')
	# ax1.plot(data[:,0], data[:,12]/data[0,12]*data[:,0]**2/data[0,0]**2, 'grey', label='68 eV')
	ax1.set(xlabel='Radius [cm]', ylabel='Np(r)*r^2 / Np(r0)*r0^2', title='Photon propagation,  max at 4th and 5th cell')
	# ax.set_ylim(1.0e-10)
	# #ax.set_xlim(3.0, 4.0)
	legend=ax1.legend(loc='center right', prop={'size': 16})
	ax1.grid()

	plt.show()
	#fig1.savefig('photonsTimesRadius2_9Bins.png')


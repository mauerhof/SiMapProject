import matplotlib.pyplot as plt
from astro import constants as c
import numpy as np
import math



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




def plotSpectrum(flux, wavelength):
	fig=plt.figure()
	plt.plot(wavelength, flux)
	#plt.yscale('log')
	#plt.xscale('log')
	#plt.ylim(1e-25, 1.1e-6)
	#plt.xlim(158,3520)
	plt.show()
	#fig.savefig('bpass100_0.002_6.3Myr.png')
	

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
	
def photoCrossSectionNu(nuBins):    #From Verner et al. (1996),  H and He
	
	epsilon0 = [0.4298, 13.61, 1.720] #eV
	sigma0 = [5.475e-14, 9.492e-16, 1.369e-14] #cm^2
	P = [2.963, 3.188, 2.963]
	ya = [32.88, 1.469, 32.88]
	yw = [0, 2.039, 0]
	y0 = [0, 0.4434, 0]
	y1 = [0, 2.136, 0]
	
	
	n = len(nuBins)
	
	xx = np.empty([3, n])
	yy = np.empty([3, n])
	sigma = np.empty([3, n])
	for element in range(3):
		for i in range(n):
			if c_cgs/nuBins[i]*1e8 > lambdaIon[element]:
				xx[element, i] = 0
				yy[element, i] = 0
				sigma[element, i] = 0
			else:
				xx[element, i] = h_cgs*nuBins[i]*erg_to_ev/epsilon0[element] - y0[element]
				yy[element, i] = math.sqrt(xx[element, i]**2 + y1[element]**2)
				sigma[element, i] = sigma0[element]*((xx[element, i] - 1)**2 + yw[element]**2) * math.pow(yy[element, i], 0.5*P[element]-5.5)/math.pow(1 + math.sqrt(yy[element,i]/ya[element]), P[element])
		
	
	return sigma

	
def Blambda(wavelength, T):   #Black body radiation, wavelength in angstroms
	
	Blambda = np.empty([wavelength.size])
	for i in range(wavelength.size):
		Blambda[i] = 2*h_cgs*c_cgs**2/wavelength[i]**5/1e-8**4/math.exp(hc_cgs/T/wavelength[i]/1e-8/kb_cgs)
	
	return Blambda
	

def meanflux(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 250, 500, 911.75) -> 3 bins
	
	sigma = photoCrossSection(wavelength)
	meanflux_array1 = flux*sigma*wavelength*1e-8/hc_cgs
	meanflux_array2 = sigma*wavelength*1e-8/hc_cgs
	
	meanflux = np.empty([len(partition)-1, 3])
	for bins in range(len(partition)-1):
		for element in range(3):
			#print wavelength[(partition[bins] < wavelength) & (wavelength < partition[bins+1])]
			meanflux[-bins-1][element] = np.trapz(meanflux_array1[element][(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])/(
			np.trapz(meanflux_array2[element][(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])]))


	return meanflux
	

def csn(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 250, 500, 911.75) -> 3 bins
	
	sigma = photoCrossSection(wavelength)
	csn_array1 = flux*sigma*wavelength*1e-8/hc_cgs
	csn_array2 = flux*wavelength*1e-8/hc_cgs
	
	csn = np.empty([len(partition)-1, 3])
	for bins in range(len(partition)-1):
		for element in range(3):
			#print wavelength[(partition[bins] < wavelength) & (wavelength < partition[bins+1])]
			csn[-bins-1][element] = np.trapz(csn_array1[element][(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])/(
			np.trapz(csn_array2[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])]))


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


	return egy*erg_to_ev #in eV
	
	
def photons(flux, wavelength, partition): #flux in erg/cm^2/s/A,  wavelength in Angstrom, partition of [0 A, 911.75 A],  example : partition = (0, 227.84, 504.26, 911.75) -> 3 bins
	
	#sigma = photoCrossSection(wavelength)
	photons_array = flux*wavelength*1e-8/hc_cgs
	
	photons = np.empty([len(partition)-1])
	for bins in range(len(partition)-1):
		#print wavelength[(partition[bins] <= wavelength) & (wavelength < partition[bins+1])]
		photons[bins] = np.trapz(photons_array[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])], wavelength[(partition[bins] <= wavelength) & (wavelength <= partition[bins+1])])


	return photons #in eV	
	
	

def setupPrint(flux, wavelength, partition, photonspersecondRamses, zero=False):
	
	photonBins = photons(flux, wavelength, partition)
	photonBinsRamses = photonBins*photonspersecondRamses/np.sum(photonBins)
	
	sigmaN = csn(flux, wavelength, partition)
	sigmaE = cse(flux, wavelength, partition)
	E = egy(flux, wavelength, partition)
	
	sep = ''
	print 'rt_n_source =',
	for i in range(len(photonBinsRamses)):
		print (sep + "{:.3e}".format(photonBinsRamses[-i-1])),
		sep = ', '
		
	print '\n'
	
	sep = ''
	for i in range(len(sigmaN[:])):
		print 'group_csn(' + str(i+1) + ',:) =',
		for j in range(3):
			if not zero:
				print sep + "{:.3e}".format(sigmaN[i][j]),
			else:
				print sep + "0.0",
			sep = ', '
		sep = ''
		print ''
		
	print '\n',
		
	sep = ''
	for i in range(len(sigmaE[:])):
		print 'group_cse(' + str(i+1) + ',:) =',
		for j in range(3):
			if not zero:
				 print sep + "{:.3e}".format(sigmaE[i][j]),
			else: 
				print sep + "0.0",
			sep = ', '
		sep = ''
		print ''
	
	print '\n',
	
	sep = ''
	print 'group_egy =',
	for i in range(len(E)):
		print sep + "{:.3e}".format(E[-i-1]),
		sep = ', '
	print '\n'
	
	
	
def CloudyPrint(flux, wavelength, partition, photonspersecondRamses): #print ionizing luminosity for Cloudy
	
	totalFlux = np.trapz(flux, wavelength)
	#print totalFlux
	photonBins = photons(flux, wavelength, partition)
	intensityCloudy = totalFlux*photonspersecondRamses/np.sum(photonBins)
	print "{:.2e}".format(intensityCloudy)
	print math.log10(intensityCloudy)
	
	
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

	ssp = readRamsesSEDs('/Users/mauerhof/Documents/seds/bpass100')
	
	#mass = 1.   #number of solar mass for the spectrum (which is in erg/s/A/Msun)
	#distance = 9.14   #distance in parsec, to get a flux from the luminosity
	#distanceFactor = 4*math.pi*(distance*pc_to_cm)**2
	
	#partition = Si()
	#partition = [0, 74, 148, 227.84, 320, 412, 504.26, 640, 775, 911.75]    #Partition of the wavelength list
	#partition = [0, 227.84, 504.26, 911.75]
	partition = [0, 911.75]
	
	print 'partition =',
	print partition
	
	
	photonsRamses = 1e48
	flux = ssp['spectra'][:, 9, 1]    #9 is for the 10th age, which here is 6.3 Myr,  1 is for the 2nd metallicity, which here is 0.002 (*solar metallicity)
	wavelength = ssp['lambdaBins']
	#flux = Blambda(ssp['lambdaBins'], 1e5)
	#plotSpectrum(flux, ssp['lambdaBins'])
	
	noIonization = False
	setupPrint(flux, ssp['lambdaBins'], partition, photonsRamses, noIonization)
	
	
	#CloudyPrint(flux, ssp['lambdaBins'], partition, photonsRamses)
	#print np.trapz(flux, wavelength)
	#print math.log10(np.trapz(flux, wavelength))
	
	IonWavelength = wavelength[wavelength<911.75]
	IonNu = np.flip(c_cgs/(IonWavelength*1e-8),0)
	IonFlux = flux[wavelength<911.75]
	
	nume = np.trapz(IonFlux*IonWavelength*1e-8*photoCrossSection(IonWavelength)/hc_cgs, IonWavelength)
	deno = 4*math.pi*np.trapz(photoCrossSectionNu(IonNu)[0]/h_cgs/IonNu, IonNu)
	print nume[0]/deno/4/math.pi*erg_to_ev
	print deno
	
	
	plotSpectrum(photoCrossSectionNu(IonNu)[0], IonNu)
	
	
	#print meanflux(flux, wavelength, [0,911.75])
	
	
	



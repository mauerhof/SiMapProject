import numpy as np
from astro import constants as c
import math

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
           spectra[:,iAge,iZ] = c.Lsun_cgs*spectrum  

   return {'ZBins':ZBins, 'ageBins':ageBins, 'lambdaBins':lambdaBins,'spectra':spectra, 'nLambda':nLambda}   #spectra are in erg/s/A/Msun


def photoCS(l, element, level):     #l in Angstroms,  level starting from 1
	
	if element == 'H':
		
		if level > 1:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
		
		E0 = np.array([0.4298]) #erg 
		cs0 = np.array([5.475e-14]) #cm^2
		P = np.array([2.963])
		ya = np.array([32.88])
		yw = np.array([0.])
		y0 = np.array([0.])
		y1 = np.array([0.])
		
		limit = c.HIon
	elif element == 'He':
		
		if level > 2:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
		
		E0 = np.array([13.61, 1.720]) #erg 
		cs0 = np.array([9.492e-16, 1.369e-14]) #cm^2
		P = np.array([3.188, 2.963])
		ya = np.array([1.469, 32.88])
		yw = np.array([2.039, 0.])
		y0 = np.array([0.4434, 0.])
		y1 = np.array([2.136, 0.])

		limit = c.HeIon
	elif element == 'Si':
		
		if level > 4:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
		
		E0 = np.array([2.317e1, 2.556e0, 1.659e-1, 12.88e0])
		cs0 = np.array([2.506e-17, 4.14e-18, 5.79e-22, 6.083e-18])
		P = np.array([3.546e0, 11.91e0, 13.36e0, 3.353e0])
		ya = np.array([20.57e0, 13.37e0, 1.474e2, 1.356e6])
		yw = np.array([2.837e-1, 1.57e0, 8.626e-1, 0e0])
		y0 = np.array([1.627e-5, 6.634e0, 9.613e1, 0e0])
		y1 = np.array([4.207e-1, 1.272e-1, 6.442e-1, 0e0])
		
		limit = c.SiIon
	else:
		print ('The element that you chose either doesn\' exist or is not implemented')
		return
		
	
	
	if c.hc_AeV/l < limit[level-1] :
		return 0e0
	else :
		xx = c.hc_AeV/l/E0[level-1] - y0[level-1]
		yy = math.sqrt(xx**2 + y1[level-1]**2)
		return cs0[level-1]*((xx - 1)**2 + yw[level-1]**2) * math.pow(yy, 0.5*P[level-1]-5.5)/math.pow(1 + math.sqrt(yy/ya[level-1]), P[level-1])

	

def csn(J, Ls, l0, l1, element, level):
	
	if element == 'H':
		if level > 1:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
	elif element == 'He':
		if level > 2:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
	elif element == 'Si':
		
		if level > 4:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
	else:
		print ('The element that you chose either doesn\' exist or is not implemented')
		return
	
	if l0 > l1:
		l0, l1 = l1, l0
	
	Ls_small = Ls[(Ls < l1) & (Ls > l0)]
	J_small = J[(Ls < l1) & (Ls > l0)]
	cs = np.empty([Ls_small.size])
	for i in range(cs.size):
		cs[i] = photoCS(Ls_small[i], element, level)
	
	return np.trapz(J_small*Ls_small*cs, Ls_small)/np.trapz(J_small*Ls_small, Ls_small)



def cse(J, Ls, l0, l1, element, level):
	
	if element == 'H':
		if level > 1:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
	elif element == 'He':
		if level > 2:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
	elif element == 'Si':
		
		if level > 4:
			print ('The ionization level that you chose either doesn\' exist or is not implemented')
			return
	else:
		print ('The element that you chose either doesn\' exist or is not implemented')
		return
	
	if l0 > l1:
		l0, l1 = l1, l0
	
	Ls_small = Ls[(Ls < l1) & (Ls > l0)]
	J_small = J[(Ls < l1) & (Ls > l0)]
	cs = np.empty([Ls_small.size])
	for i in range(cs.size):
		cs[i] = photoCS(Ls_small[i], element, level)
	
	return np.trapz(J_small*cs, Ls_small)/np.trapz(J_small, Ls_small)
	
	
	

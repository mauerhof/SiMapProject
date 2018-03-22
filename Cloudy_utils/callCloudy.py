from minirats.utils.py.cellutils import py_cell_utils as cu
from matplotlib import pyplot as plt
import numpy as np
import math
import os

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
ev_to_ryd = 0.0734989

	
	
	
	

# read cells
lmax = 6 # read all levels 

idens,ipre,ixhii,ixheii,ixheiii = 1,6,7,8,9     #Density, pressure, fraction of HII, HeII, HeIII
center = [0.5,0.5,0.5]
radius = 1.0

# RT fields
numberofBins = 7
energiesEV = np.array([10.65, 14.85, 19.77, 28.35, 37.98, 48.27, 68.25])
energiesRyd = ev_to_ryd*energiesEV #ryd
energiesErg = energiesEV/erg_to_ev

NpList = []
for i in range(numberofBins):
	NpList.append(10 + i*4)

readRT=True

ramDir = '/Users/mauerhof/Documents/Ramses/stromgren'
timestep = 17  #snapshot
ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,[idens,ipre,ixhii,ixheii,ixheiii]+NpList,center,radius,readRT)



ymin = cell_pos[:,1].min()
zmin = cell_pos[:,2].min()

X = cell_pos[:,0][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]    #Extract the list of x axis points,  for y and z at 0
#print X.max()
radius_ramses = X*3.086e21 #converts to cm

data_ramses = radius_ramses
for i in range(2, 4+numberofBins+1, 1):
	data_ramses = np.column_stack((data_ramses, cells[:,i][(cell_pos[:,1] <= ymin) & (cell_pos[:,2] <= zmin)]))

#the indices of data_ramses are :  0-radius, 1-xHII, 2-xHeII, 3-xHeIII, 4-c_r*N1, 5-c_r*N2, etc

data_ramses = data_ramses[np.argsort(data_ramses[:, 0])]  #Rearrange the data to have increasing radius


for j in range(64):
	CloudyString = "hden -1 \n" + "abundances he =-1.1 li =-40 be =-40 b =-40 c =-40 n =-40 o =-40 \n" + "continue f =-40 ne =-40 na =-40 mg =-40 \n" + "continue al =-40 si =-4.46 p =-40 s =-40 cl=-40 ar=-40 k =-40 \n" + "continue ca =-40 sc =-40 ti =-40 v =-40 cr =-40 mn =-40 fe =-40 \n" + "continue co =-40 ni =-40 cu =-40 zn =-40 \n" + "element hydrogen ionization " + "{:.3e}".format(1-data_ramses[j, 1]) + " " + "{:.3e}".format(data_ramses[j,1]) + "\n" + "element helium ionization " + "{:.3e}".format(1-data_ramses[j, 2]-data_ramses[j,3]) + " " + "{:.3e}".format(data_ramses[j,2]) + " " + "{:.3e}".format(data_ramses[j,3]) + "\n"
	for i in range(numberofBins):
		CloudyString = CloudyString + "Laser, frequency = " + "{:.3e}".format(energiesRyd[i]) + " ryd rel width 0.005 \n" + "intensity total " + "{:.3e}".format(9.79183e7*data_ramses[j,4+i]*energiesErg[i]) + " linear \n"
	CloudyString = CloudyString + "constant temp 4 \n" + "sphere \n" + "stop zone 1 \n" + "iterate to convergence max=3, error=0.1 \n" + "save averages, file=\"ion"+str(j+1)+".avr\", print last iteration \n" + "ionization, hydrogen 1 \n" + "ionization, hydrogen 2 \n" + "ionization, helium 1 \n" + "ionization, helium 2 \n" + "ionization, helium 3 \n" + "ionization, silicone 1 \n" + "ionization, silicone 2 \n" + "ionization, silicone 3 \n" + "ionization, silicone 4 \n" + "ionization, silicone 5 \n" + "ionization, silicone 6 \n" + "end of averages"
	CloudyFile=open('test.in','w')
	CloudyFile.write(CloudyString)
	CloudyFile.close()

	os.system("/Users/mauerhof/Documents/c17.00/source/cloudy.exe -r test")
	print j



# hden -1
# abundances he =-1.1 li =-40 be =-40 b =-40 c =-40 n =-40 o =-40
# continue f =-40 ne =-40 na =-40 mg =-40
# continue al =-40 si =-4.46 p =-40 s =-40 cl=-40 ar=-40 k =-40 
# continue ca =-40 sc =-40 ti =-40 v =-40 cr =-40 mn =-40 fe =-40
# continue co =-40 ni =-40 cu =-40 zn =-40
# Laser, frequency = 1.47 ryd rel width 0.005
# intensity 6.63e-5 linear
# Laser, frequency = 2.94 ryd rel width 0.005
# intensity 4.54e-5 linear
# Laser, frequency = 4.41 ryd rel width 0.005
# intensity 2.45e-7 linear
# constant temp 4
# sphere
# stop zone 1
# iterate to convergence max=3, error=0.1
# save averages, file="ion.avr", print last iteration
# ionization, hydrogen 1
# ionization, hydrogen 2
# ionization, helium 1
# ionization, helium 2
# ionization, helium 3
# end of averages

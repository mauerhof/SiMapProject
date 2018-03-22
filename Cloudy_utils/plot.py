import numpy as np
from matplotlib import pyplot as plt

kpc = 3.09e21

ion = np.genfromtxt('./ion1.avr', unpack=True)
for i in range(2,65,1):
	ion = np.column_stack((ion, np.genfromtxt('./ion'+str(i)+'.avr', unpack=True)))

radius = np.zeros_like(ion[0])
for i in range(radius.size):
	radius[i] = (i/64.+1./128.)*kpc



fig1, ax1 = plt.subplots(figsize=(9,7))
# ax1.plot(radius, ion[0], 'blue', label='HI')  
# ax1.plot(radius, ion[1], 'b--', label='HII') 
# ax1.plot(radius, ion[2], 'g', label='HeI')
# ax1.plot(radius, ion[3], 'g--', label='HeII') 
# ax1.plot(radius, ion[4], 'g:', label='HeIII') 
ax1.plot(radius, ion[5], 'r', label='SiI Cloudy') 
ax1.plot(radius, ion[6], 'r--', label='SiII Cloudy') 
ax1.plot(radius, ion[7], 'r:', label='SiIII Cloudy') 
ax1.plot(radius, ion[8], 'orange', label='SiIV Cloudy') 
ax1.plot(radius, ion[9], 'violet', label='SiV Cloudy')
ax1.set(xlabel='Radius [cm]', ylabel='Silicone fractions', title='Stromgren sphere')
# ax.set_ylim(1.0e-10)
# #ax.set_xlim(3.0, 4.0)
legend=ax1.legend(loc='center right', prop={'size': 16})
ax1.grid()

plt.show()

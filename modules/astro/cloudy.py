import math
from astro import constants as c



def writeSetup_onlySi(nH, nSi, xHII, xHeII, xHeIII, laserE, laserI, temperature, rank):
	
	CloudyString = "print off \n" + "hden " + "{:.3e}".format(math.log10(nH)) + "\n" + "elements read \n" + "helium \n" + "silicon \n" + "end of elements \n" + "element Lithium off \n" + "element Beryllium off \n" + "element Boron off \n" + "element Carbon off \n" + "element Nitrogen off \n" + "element Oxygen off \n" + "element Fluorine off \n" + "element Neon off \n" + "element Sodium off \n" + "element Magnesium off \n" + "element Aluminium off \n" + "element Phosphor off \n" + "element Sulphur off \n" + "element Chlorine off \n" + "element Argon off \n" + "element Potassium off \n" + "element Calcium off \n" + "element Scandium off \n" + "element Titanium off \n" + "element Vanadium off \n" + "element Chromium off \n" + "element Manganese off \n" + "element Iron off \n" + "element Cobalt off \n" + "element Nickel off \n" + "element Copper off \n" + "element Zinc off \n" + "abundances he =" + "{:.3e}".format(math.log10(nHe/nH)) + " si =" + "{:.3e}".format(math.log10(nSi/nH)) + " \n" + "element hydrogen ionization " + "{:.3e}".format(1.-xHII) + " " + "{:.3e}".format(xHII) + "\n" + "element helium ionization " + "{:.3e}".format(1.-xHeII - xHeIII) + " " + "{:.3e}".format(xHeII) + " " + "{:.3e}".format(xHeIII) + "\n"
	for i in range(len(laserE)):
		CloudyString = CloudyString + "Laser, frequency = " + "{:.3e}".format(laserE[i]*c.ev_to_ryd) + " ryd rel width 0.005 \n" + "intensity total " + "{:.3e}".format(laserI[i]) + " linear \n"
	CloudyString = CloudyString + "no level2 \n" + "no opacity reevaluation \n" + "no grain physics \n" + "no line transfer \n" + "database h-like levels small \n" + "database he-like levels small \n" + "no fine opacities \n" + "constant temp " + "{:.3e}".format(math.log10(temperature)) + " \n" + "sphere \n" + "stop zone 1 \n" + "c iterate to convergence max=3, error=0.1 \n" + "c print last iteration \n" + "save averages, file=\"ion" + str(rank) + ".avr\"" + "\n" + "ionization, silicone 1 \n" + "ionization, silicone 2 \n" + "ionization, silicone 3 \n" + "ionization, silicone 4 \n" + "ionization, silicone 5 \n" + "end of averages \n"# + "print on"

	CloudyFile=open("setup"+str(rank)+".in",'w')
	CloudyFile.write(CloudyString)
	CloudyFile.close()




def writeSetup(nH, Z, intensities, temperature, number=0, directory='./', printoff=False):
	
	CloudyString = ''
	if(printoff):
		CloudyString = CloudyString + 'print off \n'
	CloudyString = CloudyString + 'abundances GASS10 \n'
	CloudyString = CloudyString + 'hden ' + str(math.log10(nH)) + '\n'
	CloudyString = CloudyString + 'element abundance helium -1.103 \n'
	CloudyString = CloudyString + 'metals ' + str(math.log10(Z/0.02)) + ' log \n'				#The number is the log of the scale factor by which the abundances of all metals are multiplied
	CloudyString = CloudyString + 'element limit off '  + str(-5 + math.log10(Z/0.02)) +' \n'		#So that only  Carbon  Nitrogen  Oxygen  Neon  Magnesium  Silicon  Sulphur  Iron       
	CloudyString = CloudyString + 'table SED "sed_cloudy0.sed" \n'
	CloudyString = CloudyString + 'intensity linear ' + str(intensities[0]) + ' range 0.599 to 0.9996 Ryd \n'
	for i in range(1,len(intensities)):
		CloudyString = CloudyString + 'table SED "sed_cloudy'+str(i)+'.sed" \n'
		CloudyString = CloudyString + 'intensity linear total ' + str(intensities[i]) + ' \n'
	#CloudyString = CloudyString + 'black 1000000, luminosity = 38 \n'
	CloudyString = CloudyString + 'CMB z=0 \n'
	CloudyString = CloudyString + 'const temp ' + str(math.log10(temperature)) + '\n'
	#CloudyString = CloudyString + 'sphere \n'
	CloudyString = CloudyString + 'iterate to convergence \n'
	CloudyString = CloudyString + 'print last iteration \n'
	CloudyString = CloudyString + 'stop zone 1 \n'
	CloudyString = CloudyString + 'radius 30 \n'
	CloudyString = CloudyString + 'no molecules \n'
	CloudyString = CloudyString + 'no H2 molecule \n'
	CloudyString = CloudyString + 'no induced processes \n'
	CloudyString = CloudyString + 'save last ionization means file = "ion'+format(number,"07")+'.avr" \n'
	#CloudyString = CloudyString + 'save lines, emissivity, emergent, "lines'+format(number,"07")+'.dat" \n'
	#CloudyString = CloudyString + 'O  2 1130A \n'
	#CloudyString = CloudyString + 'end of lines \n'
	CloudyString = CloudyString + 'stop temperature off \n'
	CloudyString = CloudyString + 'set dr 0 \n'
	
	CloudyFile = open(directory+'input' + format(number,"07") + '.in', 'w')
	CloudyFile.write(CloudyString)
	CloudyFile.close()



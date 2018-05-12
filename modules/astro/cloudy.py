import math
from astro import constants as c

#BIG QUESTION :  MOLECULES ?


def writeSetup(hden, he, si, xHII, xHeII, xHeIII, laserE, laserI, temperature, rank):
	
	CloudyString = "print off hide \n" + "hden " + "{:.3e}".format(math.log10(hden)) + "\n" + "elements read \n" + "helium \n" + "silicon \n" + "end of elements \n" + "element Lithium off \n" + "element Beryllium off \n" + "element Boron off \n" + "element Carbon off \n" + "element Nitrogen off \n" + "element Oxygen off \n" + "element Fluorine off \n" + "element Neon off \n" + "element Sodium off \n" + "element Magnesium off \n" + "element Aluminium off \n" + "element Phosphor off \n" + "element Sulphur off \n" + "element Chlorine off \n" + "element Argon off \n" + "element Potassium off \n" + "element Calcium off \n" + "element Scandium off \n" + "element Titanium off \n" + "element Vanadium off \n" + "element Chromium off \n" + "element Manganese off \n" + "element Iron off \n" + "element Cobalt off \n" + "element Nickel off \n" + "element Copper off \n" + "element Zinc off \n" + "abundances he =" + "{:.3e}".format(math.log10(he)) + " si =" + "{:.3e}".format(math.log10(si)) + " \n" + "element hydrogen ionization " + "{:.3e}".format(1.-xHII) + " " + "{:.3e}".format(xHII) + "\n" + "element helium ionization " + "{:.3e}".format(1.-xHeII - xHeIII) + " " + "{:.3e}".format(xHeII) + " " + "{:.3e}".format(xHeIII) + "\n"
	for i in range(len(laserE)):
		CloudyString = CloudyString + "Laser, frequency = " + "{:.3e}".format(laserE[i]*c.ev_to_ryd) + " ryd rel width 0.005 \n" + "intensity total " + "{:.3e}".format(laserI[i]) + " linear \n"
	CloudyString = CloudyString + "no level2 \n" + "no opacity reevaluation \n" + "no grain physics \n" + "no line transfer \n" + "database h-like levels small \n" + "database he-like levels small \n" + "no fine opacities \n" + "constant temp " + "{:.3e}".format(math.log10(temperature)) + " \n" + "sphere \n" + "stop zone 1 \n" + "c iterate to convergence max=3, error=0.1 \n" + "print last iteration \n" + "save averages, file=\"ion" + str(rank) + ".avr\"" + ", print last iteration \n" + "ionization, silicone 1 \n" + "ionization, silicone 2 \n" + "ionization, silicone 3 \n" + "ionization, silicone 4 \n" + "ionization, silicone 5 \n" + "end of averages \n"# + "print on"

	CloudyFile=open("setup"+str(rank)+".in",'w')
	CloudyFile.write(CloudyString)
	CloudyFile.close()




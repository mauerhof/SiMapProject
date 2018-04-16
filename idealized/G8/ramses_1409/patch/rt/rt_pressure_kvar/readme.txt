Patch version 3:

Changed phAbs and phSc lines in rt_cooling_module.f90

In cooling_fine.f90:
	-Removed the istropization of flux in trapping
	-Making sure the trapped fraction is between 0 and 1

Version 4:
	-Trapped fraction is now exp[-2/(3*tau)]
	-Also resetting the reduced flux after trapping

Version 5:
	-Dropped the isotropic pressure, except for the trapped photons.
	Instead renormalizing the photon flux vector to cN for the pressure.

Version 6, 11 Jun 2014
	-Added dust-blackbody interaction.
	-Still need to clean up and replace rt_pressure patch

Version 7, 18 Jul 2014
	-Implicit scheme for dust/radiation temperature evolution
	-Lorentz boost for relativistic frame effects

Version 8, 1 September 2014
	-Fixed trapped photon pressure, using new nonthermal pressure variable



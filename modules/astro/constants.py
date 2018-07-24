import numpy as np

h_cgs = 6.62607e-27
c_cgs = 2.99792e10
hc_AeV = 12398.4  #A * eV
hc_cgs = h_cgs*c_cgs
kb_cgs = 1.38065e-17
Rsun_cgs = 6.957e10
ua_to_cm = 1.496e13
wien_cgs = 2.89777e-1
#numax/T is the constant wiennu*
wiennu_cgs = 5.879e10
pc_to_cm = 3.08568e18
kpc_to_cm = 1000*pc_to_cm
myr_to_s = 3.1556926e13
ryd_to_erg = 2.17987e-11
ly_to_cm = 9.461e17
me_cgs = 9.10938e-28
mp_cgs = 1.672621898e-24
epsilon0_cgs = 8.85419e-21
qe_cgs = 1.60218e-19
ev_to_erg = 1.60218e-12
A_to_cm = 10e-8
Lsun_cgs = 3.84e33
ev_to_ryd = ev_to_erg/ryd_to_erg
ryd_to_Hz = 3.28984e15

SiI, HI, SiII, HeI, SiIII, SiIV, HeII = 0,1,2,3,4,5,6


HIon = np.array([13.59843])

HeIon = np.array([24.58739, 54.417765])

SiIon = np.array([8.15168, 16.34585, 33.493, 45.14179, 166.767])

CIon = np.array([11.260288,24.3845,47.8878])

OIon = np.array([13.618055,35.12112,54.93554,77.4135,113.899])

MgIon = np.array([7.646236,15.035271,80.1436])

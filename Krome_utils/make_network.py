import numpy as np
#Data
elements_name = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P']


#Data for radiative recombination, from Badnell
Z_RR, N_RR, M_RR, W_RR, A, B, T0, T1, C, T2 = np.genfromtxt('/Users/mauerhof/Documents/rates/Radiative_rec/clist', unpack=True, skip_header=3) 


#Data for dielectronic recombination, from Badnell
Z_DR, N_DR, M_DR, W_DR = np.genfromtxt('/Users/mauerhof/Documents/rates/Dielectronic/clist_c', unpack=True, skip_header=3, usecols=(0,1,2,3))
C_DR = np.genfromtxt('/Users/mauerhof/Documents/rates/Dielectronic/clist_c', unpack=True, skip_header=3, usecols=(4, 5, 6, 7, 8, 9, 10, 11, 12))
E_DR = np.genfromtxt('/Users/mauerhof/Documents/rates/Dielectronic/clist_E', unpack=True, skip_header=3, usecols=(4, 5, 6, 7, 8, 9, 10, 11, 12))

#Data for collisional ionization, from Voronov (taken from Krome database)
dE,P,A_vor,X,K = np.genfromtxt('/Users/mauerhof/Documents/rates/Collisions/Voronov/clist', unpack=True)


def RR_Badnell(Z,N,M=1):


	if(N<0 or N>Z-1):
		print 'Error, impossible number ' + str(N) + ' of bounded electrons before recombination, for element ' + elements_name[Z-1]
		return None
	else:

		M_max = np.amax(M_RR[(Z_RR==Z)&(N_RR==N)])
		# ~ print M_max
		if M>M_max:
			M=M_max
		# ~ print 'M= ' +str(M)

		RR = str(A[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '*( sqrt(T/' + str(T0[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + ') * (1+sqrt(T/' + str(T0[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '))**('
		if(C[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]==0e0):
			RR = RR + '1-' + str(B[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0])
		else:
			RR = RR + '1-(' + str(B[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '+' + str(C[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '*exp(-' + str(T2[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '/T))'
		RR = RR + ') * (1+sqrt(T/' + str(T1[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '))**('
		if(C[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]==0e0):
			RR = RR + '1+' + str(B[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0])
		else:
			RR = RR + '1+(' + str(B[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '+' + str(C[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '*exp(-' + str(T2[(Z_RR==Z)&(M_RR==M)&(N_RR==N)][0]) + '/T))'
		RR = RR + ') )**(-1)'

		RR = RR.replace('e-','d-')
		RR = RR.replace('e+','d+')

	return RR



def DR_Badnell(Z,N,M=1):
    
	if(N<0 or N>Z-1):
		print 'Error, impossible number ' + str(N) + ' of bounded electrons before recombination, for element ' + elements_name[Z-1]
		return None
	elif(N==0):
		return '0d0'
	else:
		
		if(Z==1):
			M=1
		else:
			# ~ print M_DR[(Z_DR==Z)&(N_DR==N)]
			M_max = np.amax(M_DR[(Z_DR==Z)&(N_DR==N)])
			if M>M_max:	
				M=M_max
		
		DR = 'T**(-1.5)*(' + str(C_DR[0][(Z_DR==Z)&(M_DR==M)&(N_DR==N)][0]) + '*exp(-' + str(E_DR[0][(Z_DR==Z)&(M_DR==M)&(N_DR==N)][0]) + '/T)'
		i=1
		while((C_DR[i][(Z_DR==Z)&(M_DR==M)&(N_DR==N)][0] != 0e0) & (i<9)):
			DR = DR + ' + ' + str(C_DR[i][(Z_DR==Z)&(M_DR==M)&(N_DR==N)][0]) + '*exp(-' + str(E_DR[i][(Z_DR==Z)&(M_DR==M)&(N_DR==N)][0]) + '/T)'
			i = i+1
		DR = DR + ')'
		
		DR = DR.replace('e-','d-')
		DR = DR.replace('e+','d+')
		
		return DR
		

def Col_Voronov(Z,N):
    
	if(N<=0 or N>Z):
		print 'Error, impossible number ' + str(N) + ' of bounded electrons before collision, for element ' + elements_name[Z-1]
		return None
		
	elif(Z < 3):
		return 'auto'

	else:

		Col = str(A_vor[(Z-1)*Z/2 + Z-N - 3]) + '*(1 + ' + str(P[(Z-1)*Z/2 + Z-N - 3]) + '*sqrt(' + str(dE[(Z-1)*Z/2 + Z-N - 3]) + '*1.16045d4/T)) / (' + str(X[(Z-1)*Z/2 + Z-N - 3]) + '+' + str(dE[(Z-1)*Z/2 + Z-N - 3]) + '*1.16045d4/T) * (' + str(dE[(Z-1)*Z/2 + Z-N - 3]) + '*1.16045d4/T)**' + str(K[(Z-1)*Z/2 + Z-N - 3]) + ' * exp(-' + str(dE[(Z-1)*Z/2 + Z-N - 3]) + '*1.16045d4/T)'
		
		
		Col = Col.replace('e-','d-')
		Col = Col.replace('e+','d+')
		
		return Col



####################################################

# ~ print DR_Badnell(2,1,1)


#User section
# ~ elements_Z = range(1,15)
# ~ elements_Z = [1,2,14]
elements_Z = [1,2,8]
# ~ elements_Z = [14]
# ~ file_name = 'react_All_py'
# ~ file_name = 'react_Si_py'
file_name = 'react_O'

max_reac = 3



####################################################


counter=1
network = ''


#Recombinations
network += '@var: T = Tgas \n \n \n#RECOMBINATION, from Badnell website : \n#http://amdpp.phys.strath.ac.uk/tamoc/RR \n#http://amdpp.phys.strath.ac.uk/tamoc/DR \n \n@format:idx,R,R,P,rate \n'
for i in elements_Z:
	for j in range(min(i,max_reac)):
		network += str(counter) + ',' + elements_name[i-1] + '+'*(j+1) + ',E,' + elements_name[i-1] + '+'*j
		if(i>2):
			network += ', ' + RR_Badnell(i,i-j-1) + ' + ' + DR_Badnell(i,i-j-1) + '\n'
		else:
			# ~ network += ', ' + RR_Badnell(i,i-j-1) + ' + ' + DR_Badnell(i,i-j-1) + '\n'
			# ~ network += ', auto \n'
			network += ', 0d0 \n'
		counter += 1
	network += '\n'
	
dE,P,A_vor,X,K
	
#Collisions
network += '\n \n#COLLISIONS, from Voronov 1997 \n \n@format:idx,R,R,P,P,P,rate \n'
for i in elements_Z:
	for j in range(min(i,max_reac)):
		network += str(counter) + ',' + elements_name[i-1] + '+'*j + ',E,' + elements_name[i-1] + '+'*(j+1)  + ',E,E, '
		if(i>2):
			network += Col_Voronov(i,i-j) + '\n'
		else:
			# ~ network += 'auto \n'
			network += '0d0 \n'
		counter += 1
	network += '\n'

network += '\n \n#Photoionization \n \n@photo_start \n@format:idx,R,P,P,rate \n'
for i in elements_Z[2::]:
	for j in range(min(i,max_reac)):
		network += str(counter) + ',' + elements_name[i-1] + '+'*j + ',' + elements_name[i-1] + '+'*(j+1) + ',E, auto \n'
		counter += 1
network += '\n@photo_stop \n'


# ~ #Charge transfer for helium and silicone	
# ~ network += '\n#CHARGE TRANSFER \n \n#a*(T/1d4)**b*(1 + c*exp(d*T/1d4)) \n#Data from Kingdon & Ferland 1996 \n \n@format:idx,R,R,P,P,Tmin,Tmax,rate \n'
# ~ network += str(counter) + ',He+,H,He,H+,6d3,1d5,7.47d-15*(T/1d4)**2.06*(1 + 9.93*exp(-3.89*T/1d4)) \n' + str(counter+1) + ',He++,H,He+,H+,1d3,1d7,1d-14 \n'
# ~ network += str(counter+2) + ',Si++,H,Si+,H+,1d1,1d6,1.23d-9*(T/1d4)**0.24*(1 + 3.17*exp(4.18d-3*T/1d4)) \n'
# ~ network += str(counter+3) + ',Si+++,H,Si++,H+,1d3,3d4,4.9d-10*(T/1d4)**-8.74d-2*(1 - 0.36*exp(-0.79*T/1d4)) \n'
# ~ network += str(counter+4) + ',Si++++,H,Si+++,H+,1d3,5d4,7.58d-9*(T/1d4)**0.37*(1 + 1.06*exp(-4.09*T/1d4)) \n'
# ~ network += str(counter+5) + ',Si+,H+,Si++,H,5d3,1d6, 4.1d-10*(T/1d4)**0.24*(1 + 3.17*exp(4.18d-3*T/1d4))*exp(-3.178d4/T) \n'
# ~ counter = counter+6
	


# ~ network += '\n#Charge transfer as in Mappings \n \n#chxr = 1e-9*alpha(t4) = a*(T/1d4)**b**(1 + c*exp(d*(T/1d4))) \n#chxi = 1e-9*a*((T/1d4)**b)**(1 + c*exp(d*(T/1d4)))*exp(-e4/(T/1d4)) \n#chxi = 1e-9*a*((T/1d4)**b)*(exp(-c*(T/1d4)))*exp(-e4/(T/1d4)) \n \n@format:idx,R,R,P,P,Tmin,Tmax,rate \n'
# ~ network += str(counter) + ',He+,H,He,H+,1d1,1d5, 7.47d-15*(T/1d4)**2.06d0*(1 + 9.93d0*exp(-3.89d0*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',He++,H,He+,H+,1d3,1d7, 1d-14 \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si++,H,Si+,H+,1d0,1d5, 6.77d-9*(T/1d4)**7.36d-2*(1 - 4.3d-1*exp(-1.1d-1*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si+++,H,Si++,H+,1d3,3d4, 4.9d-10*(T/1d4)**(-8.74d-2)*(1 - 3.6d-1*exp(-7.9d-1*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si++++,H,Si+++,H+,1d3,5d4, 7.58d-9*(T/1d4)**3.7d-1*(1 + 1.06d0*exp(-4.09d0*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si+++,He,Si++,He+,1d0,1d6, 1.03d-9*(T/1d4)**6d-1*(1 - 6.1d-1*exp(-1.42d0*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si++++,He,Si+++,He+,1d0,5d5, 5.75d-10*(T/1d4)**9.3d-1*(1 + 1.33d0*exp(-2.9d-1*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si,H+,Si+,H,1d0,2d5, 9.2d-12*((T/1d4)**1.15d0)*(1 + 8d-1*exp(-2.4d-1*(T/1d4))) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si+,H+,Si++,H,1d3,1d5, 2.26d-9*((T/1d4)**7.36d-2)*(1 - 4.3d-1*exp(-1.1d-1*(T/1d4)))*exp(-3.1882d0/(T/1d4)) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si,He+,Si+,He,1d1,1d4, 1.3d-9 \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si+,He+,Si++,He,2.511d3,3d5, 1.5d-10*((T/1d4)**2.4d-1)*exp(-8.0187d0/(T/1d4)) \n'
# ~ counter += 1
# ~ network += str(counter) + ',Si++,He+,Si+++,He,3.162d3,3d5, 1.15d-9*((T/1d4)**4.4d-1)*exp(-1.0335d1/(T/1d4)) \n'
# ~ counter += 1



f=open(file_name,'w')
f.write(network)
f.close()



# ~ network += '\n@common:user_phrate \n@format:idx,R,P,P,rate \n \n'
# ~ for i in elements_Z:
	# ~ network += str(counter) + ',' + elements_name[i-1] + ',' + elements_name[i-1] + '+,E, user_phrate \n'
	# ~ counter += 1


# ~ #Collisions
# ~ network += '\n \n#COLLISIONS, auto, which corresponds to Voronov 1997 \n \n@format:idx,R,R,P,P,P,rate \n'
# ~ for i in elements_Z:
	# ~ for j in range(i):
		# ~ network += str(counter) + ',' + elements_name[i-1] + '+'*j + ',E,' + elements_name[i-1] + '+'*(j+1) 
		# ~ network += ',E,E, auto \n'
		# ~ counter += 1
	# ~ network += '\n'

